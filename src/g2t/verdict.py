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
-> likely_misidentified (reason layers_disagree: the name is wrong either way; 0 of 297 baseline
records, 23 of 23 in injections, benchmark 2026-11-21). ``strict=False`` restores the older count of method names.

Classes: likely_contaminant, likely_misidentified, numt_suspect, suspect (one family),
clean (no flag and BLAST confirms the declared name at the finest rank the record declares and the marker can
resolve), coarse_agree (BLAST agrees only at a coarser rank, e.g. family while a species is declared; see
confirmed_rank), unchecked (no positive evidence).
"""

from __future__ import annotations

import pandas as pd

COLS = ["flag", "why", "accession", "gene", "organism", "verdict", "layers", "reasons", "confirmed_rank"]
CHECK_VERDICTS = ("likely_contaminant", "likely_misidentified", "suspect", "numt_suspect")
BLAST_GATE_CHECK = {"no_hit": "BLAST 无命中", "no_close_reference": "BLAST 无足够相近的参考序列"}

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


def _flag_why(v: str, reasons: list[str], rec: dict | None, taxon: str, ran: tuple[str, ...],
              unverified: list[dict], bl: list[dict]) -> tuple[str, str]:
    """CHECK / OK / NA and the reason. CHECK means 'a person should look', not 'this is contaminated'."""
    why: list[str] = []
    if v in CHECK_VERDICTS:
        why.append(f"{v}: " + ("; ".join(reasons) if reasons else "多层证据"))
    if "blast:reverse_complement" in reasons and v not in CHECK_VERDICTS:
        why.append("方向反了（反向互补提交），序列本身需翻转")
    if rec:
        tx = str(rec.get("taxonomy", "") or "")
        if tx and "Metazoa" not in tx:
            why.append("声明的物种不是后生动物（" + "; ".join(tx.split("; ")[:2]) + "）")
        fam = str(rec.get("family", "") or "")
        if taxon and fam and fam != taxon:
            why.append(f"声明的科 {fam} 不是研究类群 {taxon}")
    if unverified:
        why.append("NCBI 分类未经核实（ncbi_unverified），不等于确认污染")
    for r in bl:
        if r.get("gate") in BLAST_GATE_CHECK and not r.get("conflict_rank"):
            why.append(BLAST_GATE_CHECK[r["gate"]] + "（可能是新物种、非目标序列或微生物，需人工看）")
            break
    if why:
        return "CHECK", " | ".join(why)
    if v in ("clean", "coarse_agree"):
        return "OK", "BLAST 与声明的名称一致" + ("" if v == "clean" else "（仅粗分类级）")
    if "blast" not in ran:
        return "NA", "未运行 BLAST，没有参考库证据；微生物检查未做"
    return "NA", "BLAST 无明确证据（unchecked）"


def combine(kmer=None, tree=None, blast=None, trans=None, organisms: dict[str, str] | None = None,
            strict: bool = True, inventory: pd.DataFrame | None = None, taxon: str = "",
            ran: tuple[str, ...] = ()) -> pd.DataFrame:
    """One row per record. ``inventory`` (accession, gene, organism, family, taxonomy: every record of the run)
    makes records without any evidence appear too; ``taxon`` is the studied family; ``ran`` names the layers that
    were run (``blast``, ``tree``), so an empty verdict can say why it is empty. Adds ``flag`` (CHECK / OK / NA)
    and ``why``, the one-line reason to look at the record."""
    k, t, b, n = (_by_acc_gene(x) for x in (kmer, tree, blast, trans))
    organisms = dict(organisms or {})
    inv: dict[tuple[str, str], dict] = {}
    if inventory is not None and len(inventory):
        for r in inventory.fillna("").to_dict("records"):
            inv[(_key(r["accession"]), str(r.get("gene", "") or ""))] = r
            organisms.setdefault(_key(r["accession"]), r.get("organism", ""))
    keys = set(k) | set(t) | set(b) | set(n) | set(inv)
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
                v = "likely_misidentified"  # both families say the name is wrong; they differ only in how far
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
        flag, why = _flag_why(v, reasons, inv.get((a, g)) or next(
            (x for (aa, _), x in inv.items() if aa == a and not g), None), taxon, ran,
            [r for r in rows_for(k, a, g) if r["flag"] == "ncbi_unverified"], bl)
        out.append(dict(flag=flag, why=why, accession=a, gene=g, organism=org, verdict=v, layers=";".join(uniq),
                        reasons="; ".join(reasons), confirmed_rank=confirmed))
    return pd.DataFrame(out, columns=COLS)
