"""One table of sequencing platform and assembler per mitogenome accession from every source (g2t-methods).

Priority, each value tagged with its source: the GenBank record itself (COMMENT / Assembly-Data / DBLINK, from
``g2t-mitostats``) > raw-read instrument fields (``g2t-sra``, genomic runs only) > the paper that reports the
record (``g2t-litmethods``, role reports / record_reference). A paper that only cites the accession is never used.
Output ``methods_final.csv``; ``*_source`` empty means no source knew.
"""
from __future__ import annotations

import argparse

import pandas as pd

from g2t import mitostats as ms

_NONE = {"", "Unknown", "No assembly info in record", ms.NOT_STATED}


def merge(stats: pd.DataFrame, sra: pd.DataFrame | None = None, lit: pd.DataFrame | None = None,
          genome: pd.DataFrame | None = None) -> pd.DataFrame:
    s = stats.fillna("")
    out = s[["accession", "organism"]].copy()
    out["platform"] = s["seq_platforms"]
    out["seq_category"] = s["seq_category"].where(s["seq_category"] != "Unknown", "")
    out["platform_source"] = s["seq_platforms"].map(lambda x: "GenBank record" if x else "")
    out["assembler"] = s.get("assemblers", "")
    out["assembler_source"] = out["assembler"].map(lambda x: "GenBank record" if x else "")
    out["assembly_group"] = s["assembler_group"]
    idx = {a: i for i, a in zip(out.index, out.accession)}

    def put(acc, col, val, src):
        i = idx.get(acc)
        if i is not None and val and out.at[i, col] in _NONE:
            out.at[i, col] = val
            out.at[i, col.split("_")[0] + "_source" if col != "platform" else "platform_source"] = src

    if sra is not None and len(sra):
        g = sra.fillna("")
        g = g[(g.get("genomic", "") == "yes")]
        for acc, grp in g.groupby("accession"):
            plats = sorted({p for p in grp["platform"] if p})
            word = {"PACBIO_SMRT": "PacBio", "OXFORD_NANOPORE": "Nanopore", "ILLUMINA": "Illumina"}
            put(acc, "platform", "; ".join(word.get(p, p) for p in plats),
                f"SRA/ENA genomic runs ({grp['run_link'].iloc[0]})")
            cat = ms.sequencing_category(" ".join(grp["platform"] + " " + grp["model"]))
            i = idx.get(acc)
            if i is not None and not out.at[i, "seq_category"] and cat:
                out.at[i, "seq_category"] = cat
    if lit is not None and len(lit):
        L = lit.fillna("")
        L = L[L.get("lit_role", "").isin(["reports", "record_reference"])]
        for _, r in L.iterrows():
            src = f"paper PMID {r.get('lit_pmid', '')} ({r.get('lit_role', '')})"
            put(r["accession"], "platform", r.get("lit_seq_platforms", ""), src)
            put(r["accession"], "assembler", r.get("lit_assemblers", ""), src)
    out["platform_source"] = out["platform_source"].where(out["platform"] != "", "")
    # nuclear genome assembly of the same BioSample: context only, never written into `assembler`
    out["genome_assembly"] = out["genome_method"] = out["genome_seq_tech"] = ""
    if genome is not None and len(genome):
        G = genome.fillna("")
        G = G[G["genome_assembly"] != ""]
        for acc, grp in G.groupby("accession"):
            i = idx.get(acc)
            if i is None:
                continue
            out.at[i, "genome_assembly"] = "; ".join(grp["genome_assembly"])
            out.at[i, "genome_method"] = "; ".join(dict.fromkeys(m for m in grp["genome_method"] if m))
            out.at[i, "genome_seq_tech"] = "; ".join(dict.fromkeys(m for m in grp["genome_seq_tech"] if m))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="g2t-methods", description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("mitostats")
    ap.add_argument("--sra", help="sra_runs.csv from g2t-sra")
    ap.add_argument("--lit", help="litmethods.csv from g2t-litmethods")
    ap.add_argument("--genome", help="genome_assemblies.csv from g2t-genome (context columns genome_*)")
    ap.add_argument("-o", "--out", default="methods_final.csv")
    a = ap.parse_args(argv)
    rd = lambda p: pd.read_csv(p, dtype=str) if p else None  # noqa: E731
    res = merge(rd(a.mitostats), rd(a.sra), rd(a.lit), rd(a.genome))
    res.to_csv(a.out, index=False)
    n = len(res)
    n_plat = (res.platform != "").sum()
    n_asm = (~res.assembler.isin(_NONE)).sum()
    print(f"{n} records: platform known {n_plat}, assembler known {n_asm} -> {a.out}")
    print(res.assembler_source.value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
