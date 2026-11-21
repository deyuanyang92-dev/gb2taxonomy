"""Unified verdict per record: combine the evidence layers (k-mer, tree, BLAST, translation).

Evidence layers (all optional DataFrames keyed by ``accession`` + ``gene``):
  kmer   screen.csv flags (near_other_genus, same_species_far, ncbi_unverified, remote)
  tree   phylo_flags.csv (outside_family, inside_other_genus, ...)
  blast  blastcheck rows (conflict_rank, note, gate)
  trans  translate_check rows (flag == numt_suspect)

Classes: likely_contaminant (>= 2 independent layers say 'another family / phylum'), likely_misidentified (>= 2 layers
say 'another genus / family' of the same phylum), numt_suspect, suspect (one layer), taxon_question (whole-genus
notes), clean (no flag and a positive BLAST / tree confirmation), unchecked (no positive evidence).
"""

from __future__ import annotations

import pandas as pd

COLS = ["accession", "gene", "organism", "verdict", "layers", "reasons"]


def _key(a: str) -> str:
    return str(a).split(".")[0]


def _by_acc(df: pd.DataFrame | None) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    if df is None or df.empty:
        return out
    for r in df.fillna("").to_dict("records"):
        out.setdefault(_key(r["accession"]), []).append(r)
    return out


def combine(kmer=None, tree=None, blast=None, trans=None, organisms: dict[str, str] | None = None) -> pd.DataFrame:
    k, t, b, n = _by_acc(kmer), _by_acc(tree), _by_acc(blast), _by_acc(trans)
    organisms = organisms or {}
    accs = sorted(set(k) | set(t) | set(b) | set(n))
    rows = []
    for a in accs:
        layers, reasons, far, near = [], [], 0, 0
        gene = ""
        org = organisms.get(a, "")
        for r in k.get(a, []):
            gene = gene or r.get("gene", "")
            org = org or r.get("organism", "")
            if r["flag"] in ("near_other_genus", "remote"):
                layers.append("kmer")
                reasons.append(f"kmer:{r['flag']}")
                far += r["flag"] == "remote"
                near += r["flag"] == "near_other_genus"
            elif r["flag"] == "same_species_far":
                layers.append("kmer")
                reasons.append("kmer:same_species_far")
                near += 1
        for r in t.get(a, []):
            gene = gene or r.get("gene", "")
            org = org or r.get("organism", "")
            if r["flag"] in ("outside_family", "inside_other_family"):
                far += 1
                layers.append("tree")
                reasons.append(f"tree:{r['flag']}")
            elif r["flag"].startswith("inside_other_genus"):
                near += 1
                layers.append("tree")
                reasons.append(f"tree:{r['flag']}")
        for r in b.get(a, []):
            gene = gene or r.get("gene", "")
            org = org or r.get("organism", "")
            cr = r.get("conflict_rank", "")
            if cr:
                layers.append("blast")
                reasons.append(f"blast:{cr}")
                if cr in ("phylum", "class", "order", "family"):
                    far += 1
                else:
                    near += 1
        numt = any(r.get("flag") == "numt_suspect" for r in n.get(a, []))
        if numt:
            layers.append("translation")
            reasons.append("translation:numt_suspect")
        uniq = list(dict.fromkeys(layers))
        positive = any(r.get("gate") == "actionable" and not r.get("conflict_rank") and r.get("lca_rank")
                       for r in b.get(a, []))
        independent = [x for x in uniq if x != "translation"]
        if far and len(independent) >= 2:
            v = "likely_contaminant"
        elif (far or near) and len(independent) >= 2:
            v = "likely_misidentified"
        elif numt and not independent:
            v = "numt_suspect"
        elif independent:
            v = "suspect"
        elif positive:
            v = "clean"
        else:
            v = "unchecked"
        rows.append(dict(accession=a, gene=gene, organism=org, verdict=v, layers=";".join(uniq),
                         reasons="; ".join(reasons)))
    return pd.DataFrame(rows, columns=COLS)
