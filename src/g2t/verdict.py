"""Unified verdict per record: combine the evidence layers (k-mer, tree, BLAST, translation).

Evidence layers (all optional DataFrames keyed by ``accession`` + ``gene``):
  kmer   screen.csv flags (near_other_genus, same_species_far, ncbi_unverified, remote)
  tree   phylo_flags.csv (outside_family, inside_other_genus, ...)
  blast  blastcheck rows (conflict_rank, note, gate)
  trans  translate_check rows (flag == numt_suspect)

One row per (accession, gene): evidence from different markers of one record is never pooled, because a COI and a
16S result are different observations, not corroboration. Evidence without a gene applies to every marker of the
accession.

Independence: the tree flag already requires a k-mer distance check, so k-mer and tree are ONE family of evidence
('intrinsic': the record against the rest of the dataset); BLAST against the reference database is the other
('reference'). Two families are needed to escalate, and they must say the same kind of thing: both 'another family /
phylum' -> likely_contaminant; both 'another genus / species' -> likely_misidentified; one far and one near
-> suspect (reason layers_disagree). ``strict=False`` restores the older count of method names.

Classes: likely_contaminant, likely_misidentified, numt_suspect, suspect (one family, or disagreeing families),
clean (no flag and BLAST confirms the declared name at the finest rank the record declares and the marker can
resolve), coarse_agree (BLAST agrees only at a coarser rank, e.g. family while a species is declared; see
confirmed_rank), unchecked (no positive evidence).
"""

from __future__ import annotations

import pandas as pd

COLS = ["accession", "gene", "organism", "verdict", "layers", "reasons", "confirmed_rank"]

FAR_RANKS = ("phylum", "class", "order", "family")
FAMILY_OF = {"kmer": "intrinsic", "tree": "intrinsic", "blast": "reference"}


def _key(a: str) -> str:
    return str(a).split(".")[0]


def _by_acc_gene(df: pd.DataFrame | None) -> dict[tuple[str, str], list[dict]]:
    out: dict[tuple[str, str], list[dict]] = {}
    if df is None or df.empty:
        return out
    for r in df.fillna("").to_dict("records"):
        out.setdefault((_key(r["accession"]), str(r.get("gene", "") or "")), []).append(r)
    return out


def _genus_of(name: str) -> str:
    """'Perinereis linea' -> 'Perinereis'; a bare genus stays; open names / empty give ''."""
    w = str(name or "").strip().split()
    return w[0] if w and w[0][:1].isupper() else ""


def _true(x) -> bool:
    return x in (True, "True", "true", 1)


class _Ev:
    """Evidence collected for one (accession, gene): method names, reasons, and the families saying far / near."""

    def __init__(self, strict: bool):
        self.strict, self.layers, self.reasons = strict, [], []
        self.far_f: set[str] = set()
        self.near_f: set[str] = set()
        self.dest: dict[str, set[str]] = {}  # evidence family -> genera the record is said to belong to
        self.reverse = False  # submitted in the reverse-complement orientation (BLAST minus-strand hit)
        self.weak = False  # a BLAST conflict resting on one reference with the declared taxon absent


def _flag(ev: _Ev, layer: str, reason: str, far: bool) -> None:
    ev.layers.append(layer)
    ev.reasons.append(f"{layer}:{reason}")
    (ev.far_f if far else ev.near_f).add(FAMILY_OF[layer] if ev.strict else layer)


def combine(kmer=None, tree=None, blast=None, trans=None, organisms: dict[str, str] | None = None,
            strict: bool = True) -> pd.DataFrame:
    k, t, b, n = (_by_acc_gene(x) for x in (kmer, tree, blast, trans))
    organisms = organisms or {}
    keys = set(k) | set(t) | set(b) | set(n)
    genes_of: dict[str, set[str]] = {}
    for a, g in keys:
        if g:
            genes_of.setdefault(a, set()).add(g)
    # evidence without a gene belongs to every marker of the accession (or to ("acc", "") when there is none)
    targets = {(a, g) for a, g in keys if g} | {(a, "") for a, g in keys if not g and a not in genes_of}

    def rows_for(d, a, g):
        rs = list(d.get((a, g), []))
        if g:
            rs += d.get((a, ""), [])
        return rs

    out = []
    for a, g in sorted(targets):
        ev = _Ev(strict)
        reasons, layers, far_f, near_f = ev.reasons, ev.layers, ev.far_f, ev.near_f

        org = organisms.get(a, "")
        for r in rows_for(k, a, g):
            org = org or r.get("organism", "")
            if r["flag"] in ("near_other_genus", "remote", "same_species_far"):
                _flag(ev, "kmer", r["flag"], r["flag"] == "remote")
                if r["flag"] != "remote" and r.get("nearest_genus"):
                    ev.dest.setdefault("intrinsic", set()).add(_genus_of(r["nearest_genus"]))
        for r in rows_for(t, a, g):
            org = org or r.get("organism", "")
            if r["flag"] in ("outside_family", "inside_other_family"):
                _flag(ev, "tree", r["flag"], True)
            elif r["flag"].startswith("inside_other_genus"):
                _flag(ev, "tree", r["flag"], False)
                gs = {_genus_of(x) for x in str(r.get("neighbours", "")).split(";")} - {""}
                if gs:
                    ev.dest.setdefault("intrinsic", set()).update(gs)
        bl = rows_for(b, a, g)
        for r in bl:
            org = org or r.get("organism", "")
            cr = r.get("conflict_rank", "")
            if cr and _true(r.get("weak")):
                ev.weak = True
                reasons.append(f"blast:{cr} (weak: one reference, declared taxon absent from hits)")
            elif cr:
                _flag(ev, "blast", cr, cr in FAR_RANKS)
                if cr in ("genus", "species") and _genus_of(r.get("lca_taxon", "")):
                    ev.dest.setdefault("reference", set()).add(_genus_of(r["lca_taxon"]))
            if _true(r.get("reverse")):
                reasons.append("blast:reverse_complement")
                ev.reverse = True
        numt = any(r.get("flag") == "numt_suspect" for r in rows_for(n, a, g))
        if numt:
            reasons.append("translation:numt_suspect")
        uniq = list(dict.fromkeys(layers))
        fams = far_f | near_f
        agreeing = [r for r in bl if r.get("gate") == "actionable" and not r.get("conflict_rank")
                    and _true(r.get("agrees"))]
        # older blastcheck files have no 'resolved' column: agreement then counts as confirmation
        resolved = any(_true(r.get("resolved", r.get("agrees"))) for r in agreeing)
        confirmed = next((r.get("confirmed_rank", "") for r in agreeing if r.get("confirmed_rank")), "")
        if len(fams) >= 2:
            if not near_f:
                v = "likely_contaminant"
            elif not far_f:
                v = "likely_misidentified"
            else:
                v = "suspect"
                reasons.append("layers_disagree: one says another family, the other another genus / species")
        elif fams or ev.weak:
            v = "suspect"
        elif numt:
            v = "numt_suspect"
        elif resolved:
            v = "clean"
        elif agreeing:
            v = "coarse_agree"
        else:
            v = "unchecked"
        if ev.reverse and v in ("clean", "coarse_agree", "unchecked"):
            v = "suspect"  # wrong orientation is a record error even when the name checks out
        d = ev.dest
        if strict and v == "likely_misidentified" and d.get("intrinsic") and d.get("reference") \
                and not d["intrinsic"] & d["reference"]:
            # two families agree that the name is wrong but not on where the record belongs
            v = "suspect"
            reasons.append("destinations_disagree: " + "/".join(sorted(d["intrinsic"])) + " vs "
                           + "/".join(sorted(d["reference"])))
        if not strict and len(fams) >= 2:  # old rule: any two layers, far and near pooled
            v = "likely_contaminant" if far_f else "likely_misidentified"
        out.append(dict(accession=a, gene=g, organism=org, verdict=v, layers=";".join(uniq),
                        reasons="; ".join(reasons), confirmed_rank=confirmed))
    return pd.DataFrame(out, columns=COLS)
