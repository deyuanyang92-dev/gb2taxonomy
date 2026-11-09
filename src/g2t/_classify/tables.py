"""Record filters, filtered/status file names, output files and column order."""

from __future__ import annotations

import contextlib
import logging
import os

import pandas as pd

from g2t.utils import parse_interval

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


def filter_dataframe(
    df: pd.DataFrame, moleculetype: str, mol_type: str, length_filter: str, organelle_filter: str
) -> pd.DataFrame:
    if moleculetype.lower() != "all":
        df = df[df["MoleculeType"].str.lower() == moleculetype.lower()]
    if mol_type.lower() != "all":
        allowed = [x.strip().lower() for x in mol_type.split(",")]
        df = df[df["mol_type"].str.lower().isin(allowed)]
    if length_filter.lower() != "all":
        lo, hi = parse_interval(length_filter)
        if lo is not None:
            df = df[df["Length"] >= lo]
        if hi is not None:
            df = df[df["Length"] <= hi]
    if organelle_filter.lower() != "all":
        df = df[df["organelle"].str.lower() == organelle_filter.lower()]
    return df


FILTERED_FILE = "filtered_records.csv"
STATUS_FILE = "record_status.csv"
_FILTER_COLS = ["input_row", "ACCESSION", "LocusID", "Organism", "Length", "Definition", "step", "reason"]


def _write_filtered(output_dir: str, rows: list[dict], append: bool) -> None:
    """Write records removed by a filter, with the reason, so nothing disappears silently."""
    path = os.path.join(output_dir, FILTERED_FILE)
    df = pd.DataFrame(rows, columns=_FILTER_COLS)
    if append and os.path.exists(path):
        df = pd.concat([pd.read_csv(path, dtype=str), df.astype(str)], ignore_index=True)
    df.to_csv(path, index=False)


def _filtered_rows(removed: pd.DataFrame, step: str, reason) -> list[dict]:
    out = []
    for _, r in removed.iterrows():
        out.append(
            {
                "input_row": r.get("g2t_input_row", ""),
                "ACCESSION": r.get("ACCESSION", ""),
                "LocusID": r.get("LocusID", ""),
                "Organism": r.get("Organism", r.get("organism", "")),
                "Length": r.get("g2t_length_raw", "") if pd.isna(r.get("Length")) else int(r.get("Length")),
                "Definition": r.get("Definition", ""),
                "step": step,
                "reason": reason(r) if callable(reason) else reason,
            }
        )
    return out


def _filter_dataframe_logged(
    df: pd.DataFrame, moleculetype: str, mol_type: str, length_filter: str, organelle_filter: str, removed: list[dict]
) -> pd.DataFrame:
    """filter_dataframe, but every removed row is appended to `removed` with its reason."""

    def apply(d, keep, reason):
        removed.extend(_filtered_rows(d[~keep], "classify round 1", reason))
        return d[keep]

    if moleculetype.lower() != "all":
        df = apply(
            df,
            df["MoleculeType"].str.lower() == moleculetype.lower(),
            lambda r: f"MoleculeType '{r.get('MoleculeType', '')}' is not '{moleculetype}'",
        )
    if mol_type.lower() != "all":
        allowed = [x.strip().lower() for x in mol_type.split(",")]
        df = apply(
            df,
            df["mol_type"].str.lower().isin(allowed),
            lambda r: f"mol_type '{r.get('mol_type', '')}' not in {mol_type}",
        )
    if length_filter.lower() != "all":
        lo, hi = parse_interval(length_filter)
        keep = pd.Series(True, index=df.index)
        if lo is not None:
            keep &= df["Length"] >= lo
        if hi is not None:
            keep &= df["Length"] <= hi
        df = apply(df, keep, lambda r: f"length {int(r['Length'])} bp outside --length {length_filter}")
    if organelle_filter.lower() != "all":
        df = apply(
            df,
            df["organelle"].str.lower() == organelle_filter.lower(),
            lambda r: f"organelle '{r.get('organelle', '')}' is not '{organelle_filter}'",
        )
    return df


CLASSIFY_OUTPUTS = [
    "assigned_genes_type.csv",
    "assigned_genes_types2.csv",
    "assigned_genes_types_all.csv",
    "unmatched_sequences.csv",
    "unmatched_sequences2.csv",
    "unmatched_categorized.csv",
    "unmatched_report.txt",
    "conflicted_genes.csv",
    "multiple_matched.csv",
    "multiple_matched2.csv",
    FILTERED_FILE,
    STATUS_FILE,
]


def _clean_outputs(output_dir: str, extra: list[str] | None = None) -> None:
    """Remove the products of a previous run, so stale files are never read back."""
    for name in CLASSIFY_OUTPUTS + list(extra or []):
        path = os.path.join(output_dir, name)
        if os.path.isfile(path):
            os.remove(path)


def _write_all(output_dir: str, prematch_output: str) -> str:
    """assigned_genes_types_all.csv = round 1 (+ round 2 if it ran); always written."""
    frames = []
    for name in (prematch_output, "assigned_genes_types2.csv"):
        path = os.path.join(output_dir, name)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            with contextlib.suppress(pd.errors.EmptyDataError):
                frames.append(pd.read_csv(path, dtype=str))
    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=DIAG_COLUMNS + ["LocusID"])
    if "LocusID" in merged.columns:
        merged = merged.drop_duplicates(subset=["LocusID"])
    merged = reorder_columns(merged)
    all_file = os.path.join(output_dir, "assigned_genes_types_all.csv")
    merged.to_csv(all_file, index=False)
    return all_file


DIAG_COLUMNS = ["gene_type", "Conflict", "match_source", "Original_match", "Conflict_reason", "Assignment_reason"]


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    front = [c for c in DIAG_COLUMNS if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest]
