#!/usr/bin/env python3
"""Sensitivity of specimen linking to each kind of evidence: re-run steps 3b and 4 with one evidence kind
disregarded (ReconcileConfig.ignore_evidence) and score against the published tables.

    python benchmark/sensitivity.py OUT.tsv GOLD.tsv:DATASET:RUN_DIR [...]
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import evaluate as ev  # noqa: E402
import pandas as pd  # noqa: E402

from g2t.organize import organize  # noqa: E402
from g2t.reconcile import ReconcileConfig, reconcile  # noqa: E402

SCENARIOS = {
    "all evidence (default)": (),
    "without shared publication": ("same publication",),
    "without same collection date": ("same collection date",),
    "without coordinates": ("same coordinates", "coordinates"),
    "without identical identifier": ("identical voucher",),
    "without author, collector, locality": ("same first author", "same collector", "same locality"),
}


def main() -> None:
    out, specs = sys.argv[1], sys.argv[2:]
    rows = []
    for spec in specs:
        gold_path, dataset, run = spec.split(":")
        gold = pd.read_csv(gold_path, sep="\t", dtype={"specimen_level": bool})
        gold = gold[gold.dataset == dataset]
        in_matrix = set(ev.matrix_clusters(run))
        gold = gold[gold.accession.isin(in_matrix)]
        voucher_csv = Path(run) / "updated_species_vouchers" / "updated_species_voucher.csv"
        for name, ignore in SCENARIOS.items():
            with tempfile.TemporaryDirectory() as tmp:
                reconcile(str(voucher_csv), tmp, ReconcileConfig(ignore_evidence=ignore))
                Path(tmp, "organized_genes").mkdir()
                organize(
                    str(Path(tmp, "reconciled_species_voucher.csv")),
                    str(Path(tmp, "organized_genes", "organized_species_voucher.csv")),
                )
                s = ev.linking_scores(gold, ev.matrix_clusters(tmp))
            rows.append(
                {
                    "dataset": dataset,
                    "scenario": name,
                    "precision": s["precision"],
                    "recall": s["recall"],
                    "F1": s["F1"],
                    "false_pairs": s["false_merges"],
                    "specimens_recovered": f"{s['specimens_complete']}/{s['specimens']}",
                }
            )
            print(json.dumps(rows[-1]), flush=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
