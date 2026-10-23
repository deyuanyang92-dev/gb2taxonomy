#!/usr/bin/env python3
"""Score g2t output against the gold standard built by build_gold.py.

Two checks per dataset:

1. Gene type: the gene column of each accession in the published table vs. the matrix column(s) where
   g2t puts it; accessions missing from the matrix get their status from labeled_genes/record_status.csv.
2. Specimen linking (specimen-level rows only): pairwise precision/recall/F1 over accession pairs
   ("same specimen" in the published table vs. "same matrix row" in g2t), plus the share of
   published specimens whose accessions all end up in one row and no other gold accession joins them.

Linking methods compared on the same records:
    singleton   every record its own row (no linking)
    exact_id    same organism + identical identifier string (specimen_voucher, else strain, else isolate),
                i.e. linking by exact voucher labels
    supercrunch SuperCRUNCH's own parser (Parse_Loci.parse_fasta_record, Portik & Wiens 2020) applied to the
                FASTA title of each record: same species name + same voucher/isolate/strain code from the
                DEFINITION line (its "vouchered" mode); records without a code stay unlinked
                (needs --supercrunch path/to/Parse_Loci.py)
    g2t_exact   g2t matrix built with --skip_reconcile
    g2t         g2t matrix (default: voucher variants reconciled)

Usage:
    python benchmark/evaluate.py --gold benchmark/gold/gold_standard.tsv \
        --dataset Haplotaxidae --run out_haplotaxis --run-exact out_haplotaxis_norec
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

GENE_COLS_PLANT = [
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
GENE_COLS = [
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


def _base(acc) -> str:
    return str(acc).split(".")[0].strip()


def matrix_cells(run_dir: str) -> list[tuple[str, int, str]]:
    """(accession, matrix row index, gene column) for every accession in the matrix."""
    m = pd.read_csv(f"{run_dir}/organized_genes/organized_species_voucher.csv", dtype=str)
    out = []
    for i, row in m.iterrows():
        for col in GENE_COLS + GENE_COLS_PLANT:
            if col in m.columns and pd.notna(row[col]):
                out += [(_base(tok), i, col) for tok in re.split(r"[;,\s]+", row[col]) if tok]
    return out


def matrix_clusters(run_dir: str) -> dict[str, int]:
    """accession -> matrix row index."""
    return {acc: i for acc, i, _ in matrix_cells(run_dir)}


def exact_id_clusters(run_dir: str) -> dict[str, object]:
    meta = pd.read_csv(f"{run_dir}/gb_metadata/final.csv", dtype=str).fillna("")
    out: dict[str, object] = {}
    for _, r in meta.iterrows():
        ident = next((r[c].strip() for c in ("specimen_voucher", "strain", "isolate") if r.get(c, "").strip()), "")
        acc = _base(r["ACCESSION"])
        out[acc] = (r["Organism"], ident) if ident else ("single", acc)
    return out


def supercrunch_clusters(run_dir: str, parse_loci: str) -> dict[str, object]:
    import importlib.util

    spec = importlib.util.spec_from_file_location("sc_parse_loci", parse_loci)
    mod = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    meta = pd.read_csv(f"{run_dir}/gb_metadata/final.csv", dtype=str).fillna("")
    out: dict[str, object] = {}
    for _, r in meta.iterrows():
        title = f">{r.get('Version') or r['ACCESSION']} {r['Definition']}"
        _, sp, _, _, voucher = mod.parse_fasta_record(title, set(), set(), True)
        acc = _base(r["ACCESSION"])
        out[acc] = (sp, voucher) if voucher != "NA" else ("single", acc)
    return out


def exact_any_field_clusters(run_dir: str) -> dict[str, object]:
    """Same organism + an identical identifier in any voucher field (case/punctuation ignored), transitively."""
    meta = pd.read_csv(f"{run_dir}/gb_metadata/final.csv", dtype=str).fillna("")
    parent: dict[str, str] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for _, r in meta.iterrows():
        acc = _base(r["ACCESSION"])
        find(acc)
        for c in ("specimen_voucher", "isolate", "culture_collection", "clone", "strain"):
            v = re.sub(r"[^A-Z0-9]", "", str(r.get(c, "")).upper())
            if len(v) >= 3:
                parent[find(acc)] = find(f"{r['Organism']}|{v}")
    return {a: find(a) for a in list(parent) if "|" not in a}


def core_only_clusters(run_dir: str) -> dict[str, object]:
    """Same organism + any shared voucher core (g2t's candidate rule) merged WITHOUT evidence (ablation)."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from g2t.reconcile import VOUCHER_FIELDS, voucher_cores

    meta = pd.read_csv(f"{run_dir}/gb_metadata/final.csv", dtype=str).fillna("")
    parent: dict[str, str] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for _, r in meta.iterrows():
        acc = _base(r["ACCESSION"])
        find(acc)
        for c in VOUCHER_FIELDS:
            v = str(r.get(c, ""))
            if any(ch.isalnum() for ch in v):
                for core in voucher_cores(v):
                    parent[find(acc)] = find(f"{r['Organism']}|{core}")
    return {a: find(a) for a in list(parent) if "|" not in a}


def pairs(groups: dict[str, object], accs: list[str]) -> set[tuple[str, str]]:
    by = defaultdict(list)
    for a in accs:
        by[groups.get(a, ("single", a))].append(a)
    return {tuple(sorted(p)) for members in by.values() for p in itertools.combinations(members, 2)}


def specimen_gold(gold: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Single-specimen rows; accessions listed under two specimens (a table error, e.g. PP968595) are dropped."""
    spec = gold[gold.specimen_level].drop_duplicates(["row_id", "accession"])
    n_rows = spec.groupby("accession").row_id.nunique()
    ambiguous = set(n_rows[n_rows > 1].index)
    return spec[~spec.accession.isin(ambiguous)], len(ambiguous)


def bcubed(truth: dict[str, object], groups: dict[str, object], accs: list[str]) -> tuple[float, float, float]:
    """B-cubed precision/recall/F1 (Bagga & Baldwin 1998) over the gold accessions."""
    pred = {a: groups.get(a, ("single", a)) for a in accs}
    by_t, by_p = defaultdict(set), defaultdict(set)
    for a in accs:
        by_t[truth[a]].add(a)
        by_p[pred[a]].add(a)
    p = sum(len(by_t[truth[a]] & by_p[pred[a]]) / len(by_p[pred[a]]) for a in accs) / len(accs)
    r = sum(len(by_t[truth[a]] & by_p[pred[a]]) / len(by_t[truth[a]]) for a in accs) / len(accs)
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)


def bootstrap(gold: pd.DataFrame, methods: dict[str, dict], n: int = 1000, seed: int = 1) -> dict:
    """95% intervals over specimens (resampled with replacement) for pairwise F1 and B-cubed F1 of each method,
    and for the difference g2t - exact_id."""
    import random

    spec, _ = specimen_gold(gold)
    rows = sorted(spec.row_id.unique())
    by_row = {r: list(g.accession) for r, g in spec.groupby("row_id")}
    rng = random.Random(seed)
    draws: dict[str, list[float]] = defaultdict(list)
    for _ in range(n):
        pick = [rng.choice(rows) for _ in rows]
        truth, accs = {}, []
        for k, r in enumerate(pick):  # duplicated specimens get their own copies
            for a in by_row[r]:
                truth[f"{a}#{k}"] = f"{r}#{k}"
                accs.append(f"{a}#{k}")
        row_of = {x: r for k, r in enumerate(pick) for x in (f"{a}#{k}" for a in by_row[r])}
        for name, cl in methods.items():
            g = {f"{a}#{k}": cl.get(a, ("single", a)) for k, r in enumerate(pick) for a in by_row[r]}
            gp, pp = pairs(truth, accs), pairs(g, accs)
            # pairs between two copies of one resampled specimen are artefacts of resampling
            pp = {(x, y) for x, y in pp if not (row_of[x] == row_of[y] and x.split("#")[1] != y.split("#")[1])}
            tp = len(gp & pp)
            pr = tp / len(pp) if pp else 1.0
            rc = tp / len(gp) if gp else 1.0
            draws[f"{name}_F1"].append(2 * pr * rc / (pr + rc) if pr + rc else 0.0)
            # B-cubed with the same rule: an item's predicted cluster excludes copies of its own specimen
            by_pred = defaultdict(list)
            for x in accs:
                by_pred[g[x]].append(x)
            bp = br = 0.0
            for x in accs:
                cl = [y for y in by_pred[g[x]] if not (row_of[x] == row_of[y] and x.split("#")[1] != y.split("#")[1])]
                same = sum(1 for y in cl if truth[y] == truth[x])
                bp += same / len(cl)
                br += same / len(by_row[row_of[x]])
            bp, br = bp / len(accs), br / len(accs)
            draws[f"{name}_bcubed_F1"].append(2 * bp * br / (bp + br) if bp + br else 0.0)
    if "g2t_F1" in draws and "exact_id_F1" in draws:
        draws["g2t_minus_exact_id_F1"] = [a - b for a, b in zip(draws["g2t_F1"], draws["exact_id_F1"])]
        draws["g2t_minus_exact_id_bcubed_F1"] = [
            a - b for a, b in zip(draws["g2t_bcubed_F1"], draws["exact_id_bcubed_F1"])
        ]
    out = {}
    for k, v in draws.items():
        v = sorted(v)
        out[k] = [round(v[int(0.025 * n)], 4), round(v[int(0.975 * n) - 1], 4)]
    return out


def linking_scores(gold: pd.DataFrame, groups: dict[str, object]) -> dict:
    spec, n_ambiguous = specimen_gold(gold)
    accs = sorted(set(spec.accession))
    truth = {a: r for r, a in zip(spec.row_id, spec.accession)}
    gold_pairs = pairs(truth, accs)
    pred_pairs = pairs(groups, accs)
    tp = len(gold_pairs & pred_pairs)
    p = tp / len(pred_pairs) if pred_pairs else 1.0
    r = tp / len(gold_pairs) if gold_pairs else 1.0
    f = 2 * p * r / (p + r) if p + r else 0.0
    complete = 0
    for _, g in spec.groupby("row_id"):
        keys = {groups.get(a, ("single", a)) for a in g.accession}
        if len(keys) == 1:
            key = keys.pop()
            joined = {a for a in accs if groups.get(a, ("single", a)) == key}
            complete += joined == set(g.accession)
    n_spec = spec.row_id.nunique()
    return {
        "specimens": n_spec,
        "accessions": len(accs),
        "gold_pairs": len(gold_pairs),
        "pred_pairs": len(pred_pairs),
        "false_merges": len(pred_pairs - gold_pairs),
        "precision": round(p, 4),
        "recall": round(r, 4),
        "F1": round(f, 4),
        "specimens_complete": complete,
        "bcubed_precision": round(bcubed(truth, groups, accs)[0], 4),
        "bcubed_recall": round(bcubed(truth, groups, accs)[1], 4),
        "bcubed_F1": round(bcubed(truth, groups, accs)[2], 4),
        "ambiguous_accessions_dropped": n_ambiguous,
        "specimens_complete_pct": round(100 * complete / n_spec, 1) if n_spec else None,
    }


def gene_scores(gold: pd.DataFrame, run_dir: str) -> tuple[dict, pd.DataFrame]:
    st = pd.read_csv(f"{run_dir}/labeled_genes/record_status.csv", dtype=str)
    st["acc"] = st.ACCESSION.map(_base)
    st = st.drop_duplicates("acc").set_index("acc")
    g = gold.drop_duplicates(["accession", "gene"]).copy()
    cols: dict[str, set[str]] = defaultdict(set)
    for acc, _, col in matrix_cells(run_dir):
        cols[acc].add(col)

    def judge(r):
        if r.accession not in st.index:
            return "not_in_input"
        if r.accession in cols:
            return "correct" if r.gene in cols[r.accession] else "wrong"
        s = st.loc[r.accession]
        return s.status if s.status != "assigned" else "not_in_matrix"

    g["result"] = g.apply(judge, axis=1)
    g["g2t_columns"] = [",".join(sorted(cols.get(a, ()))) for a in g.accession]
    counts = g.result.value_counts().to_dict()
    inp = g[g.result != "not_in_input"]
    counts["accuracy_of_input"] = round((inp.result == "correct").mean(), 4) if len(inp) else None
    return counts, g


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--run", required=True, help="g2t output directory (default settings)")
    ap.add_argument("--run-exact", help="g2t output directory built with --skip_reconcile")
    ap.add_argument("--run-high", help="g2t output directory built with --reconcile_min_confidence high")
    ap.add_argument("--bootstrap", type=int, default=0, help="bootstrap resamples for 95%% intervals (e.g. 1000)")
    ap.add_argument("--details", help="write per-accession gene results here (TSV)")
    ap.add_argument("--supercrunch", help="path to SuperCRUNCH's Parse_Loci.py (adds the 'supercrunch' method)")
    a = ap.parse_args()

    gold = pd.read_csv(a.gold, sep="\t", dtype={"specimen_level": bool})
    gold = gold[gold.dataset == a.dataset]
    if gold.empty:  # dataset without a published table: report nothing to score
        print(json.dumps({"dataset": a.dataset, "gene_type": {}, "linking": {}}, indent=1))
        return
    genes, detail = gene_scores(gold, a.run)
    if a.details:
        detail.to_csv(a.details, sep="\t", index=False)
    result = {"dataset": a.dataset, "gene_type": genes, "linking": {}}
    if gold.specimen_level.any():
        methods = {
            "singleton": {},
            "exact_id": exact_id_clusters(a.run),
            "exact_any_field": exact_any_field_clusters(a.run),
            "core_only": core_only_clusters(a.run),
            "g2t": matrix_clusters(a.run),
        }
        if a.run_high:
            methods["g2t_high"] = matrix_clusters(a.run_high)
        if a.supercrunch:
            methods["supercrunch"] = supercrunch_clusters(a.run, a.supercrunch)
        if a.run_exact:
            methods["g2t_exact"] = matrix_clusters(a.run_exact)
        result["linking"] = {k: linking_scores(gold, v) for k, v in methods.items()}
        if a.bootstrap:
            result["bootstrap_95"] = bootstrap(
                gold, {k: v for k, v in methods.items() if k != "singleton"}, a.bootstrap
            )
    print(json.dumps(result, indent=1, default=str))


if __name__ == "__main__":
    main()
