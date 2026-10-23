#!/usr/bin/env python3
"""Collect the benchmark results into two tables (TSV + Markdown) in WORKDIR/results.

datasets.tsv  per dataset: records, status counts, specimens and multi-gene specimens (with and without
              voucher reconciliation), gene-type accuracy, run time
linking.tsv   per dataset and linking method: pairwise precision/recall/F1, false merges, complete specimens

Usage:
    python benchmark/summarize.py WORKDIR
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

GENES_EXTRA = [
    "nd1",
    "nd2",
    "nd3",
    "nd4",
    "nd4l",
    "nd5",
    "nd6",
    "atp6",
    "atp8",
    "plastome",
    "rbcl",
    "matk",
    "trnl-trnf",
    "psba-trnh",
    "rpl16",
    "rps16",
    "ndhf",
    "rpoc1",
    "ycf1",
    "rpl32-trnl",
    "trns-trng",
    "atpf-atph",
    "petb-petd",
    "psbe-petl",
    "rbcl-atpb",
    "trnk",
]
GENES = [
    "mtgenome",
    "coi",
    "16s",
    "12s",
    "cob",
    "cox2",
    "cox3",
    "nd1",
    "nd2",
    "nd3",
    "nd4",
    "nd4l",
    "nd5",
    "nd6",
    "atp6",
    "atp8",
    "18s",
    "28s",
    "its1-its2",
    "18-28s",
    "ef-1",
    "h3",
]
NAMES = {
    "haplotaxidae": "Haplotaxidae",
    "hesionidae": "Hesionidae",
    "opheliidae": "Opheliidae",
    "phyllodocida": "Phyllodocida",
    "priapulidae": "Priapulidae",
    "nymphaeaceae": "Nymphaeaceae",
    "ginkgoaceae": "Ginkgoaceae",
    "ophioglossaceae": "Ophioglossaceae",
}


def matrix_stats(run: Path) -> tuple[int, int]:
    m = pd.read_csv(run / "organized_genes" / "organized_species_voucher.csv", dtype=str)
    n_genes = m[[g for g in dict.fromkeys(GENES + GENES_EXTRA) if g in m.columns]].notna().sum(axis=1)
    return len(m), int((n_genes >= 2).sum())


def main(work: str) -> None:
    w = Path(work)
    res = w / "results"
    timing = pd.read_csv(res / "timing.tsv", sep="\t")
    rows, link = [], []
    for name, label in NAMES.items():
        f = res / f"{name}.json"
        if not f.exists():
            continue
        r = json.loads(f.read_text())
        st = pd.read_csv(w / f"out_{name}" / "labeled_genes" / "record_status.csv")
        spec, multi = matrix_stats(w / f"out_{name}")
        spec_x, multi_x = matrix_stats(w / f"out_{name}_exact")
        t = timing[timing.dataset == name].set_index("mode")
        gt = r["gene_type"]
        rows.append(
            {
                "dataset": label,
                "gold": r["dataset"],
                "records": len(st),
                "assigned": int((st.status == "assigned").sum()),
                "unmatched": int((st.status == "unmatched").sum()),
                "filtered": int((st.status == "filtered").sum()),
                "specimens_exact": spec_x,
                "multigene_exact": multi_x,
                "specimens_g2t": spec,
                "multigene_g2t": multi,
                "gold_accessions_in_input": sum(
                    v for k, v in gt.items() if k not in ("not_in_input", "accuracy_of_input")
                ),
                "gene_type_correct": gt.get("correct", 0),
                "gene_type_accuracy": gt.get("accuracy_of_input"),
                "download_s": t.seconds.get("download"),
                "g2t_s": t.seconds.get("g2t"),
                "g2t_max_rss_mb": t.max_rss_mb.get("g2t"),
            }
        )
        for method, v in r["linking"].items():
            link.append({"dataset": label, "gold": r["dataset"], "method": method, **v})
    ds, lk = pd.DataFrame(rows), pd.DataFrame(link)
    ds.to_csv(res / "datasets.tsv", sep="\t", index=False)
    lk.to_csv(res / "linking.tsv", sep="\t", index=False)
    md = ["## Datasets", "", ds.to_markdown(index=False), "", "## Specimen linking", "", lk.to_markdown(index=False)]
    (res / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main(sys.argv[1])
