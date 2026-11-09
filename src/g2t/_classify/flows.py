"""Run round 1 and round 2 over a table."""

from __future__ import annotations

import logging
import os
from typing import Any

import pandas as pd

from g2t._classify.config import MatchConfig
from g2t._classify.conflict import conflict_check_and_reassign
from g2t._classify.match import _mtgenome_min, _record_group, first_round_match
from g2t._classify.recheck import recheck_match_row
from g2t._classify.report import _save_csv
from g2t._classify.tables import (
    DIAG_COLUMNS,
    FILTERED_FILE,
    _filter_dataframe_logged,
    _filtered_rows,
    _write_filtered,
    reorder_columns,
)
from g2t._classify.text import clean_definition, clean_length
from g2t.utils import UnsupportedFormatError, parse_interval
from g2t.utils import read_table as utils_read_table

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


# =========================
# Row processing
# =========================
def process_row(row: dict[str, Any], config: MatchConfig) -> dict[str, Any]:
    prelim = first_round_match(row, config)
    gene_type_str = prelim["gene_type"]

    if not gene_type_str:
        return prelim

    definition = clean_definition(row.get("Definition", ""))
    group = _record_group(str(row.get("organelle", "")).lower(), definition)
    length_val = clean_length(row.get("Length", "0"))

    final = conflict_check_and_reassign(
        gene_type_str,
        definition,
        "",
        group,
        length_val,
        mt_min=_mtgenome_min(config),
        cp_min=parse_interval(config.get_gene_length_range("plastome"))[0],
    )

    if final["match_source"] == "reassign_mode":
        match_source = "first_match -> conflict_check"
    else:
        match_source = prelim["match_source"]

    reasons = []
    if prelim["Assignment_reason"]:
        reasons.append(prelim["Assignment_reason"])
    if final["Assignment_reason"]:
        reasons.append(f"[conflict check] {final['Assignment_reason']}")

    return {
        "gene_type": final["gene_type"],
        "match_source": match_source,
        "Original_match": final["Original_match"],
        "conflict": final["conflict"],
        "Conflict_reason": final["Conflict_reason"],
        "Assignment_reason": "; ".join(reasons),
    }


# =========================
# Main flows
# =========================
def process_prematch(
    input_file: str,
    output_dir: str,
    config: MatchConfig,
    prematch_output: str,
    multip_output: str,
    unmatched_output: str,
    conflicted_output: str,
    moleculetype: str,
    mol_type: str,
    length_filter: str,
    organelle_filter: str,
) -> pd.DataFrame:
    os.makedirs(output_dir, exist_ok=True)

    ext = os.path.splitext(input_file)[1].lower()
    if ext not in (".csv", ".tsv", ".xls", ".xlsx", ".txt"):
        logger.error(f"Unsupported format: {ext}")
        raise UnsupportedFormatError(f"Unsupported format: {ext}")
    df = utils_read_table(input_file)
    df.columns = df.columns.str.strip()

    df["g2t_length_raw"] = df["Length"].astype(str)
    df["Length"] = df["Length"].astype(str).str.replace(r"\s*bp\s*", "", regex=True).str.replace(",", "")
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
    df["g2t_input_row"] = range(len(df))
    removed: list[dict] = []

    missing = df["Length"].isna()
    removed.extend(_filtered_rows(df[missing], "classify round 1", "length missing or not numeric"))
    df = df[~missing]
    lo, hi = parse_interval(config.global_length_range)
    keep = pd.Series(True, index=df.index)
    if lo is not None:
        keep &= df["Length"] >= lo
    if hi is not None:
        keep &= df["Length"] <= hi
    # plastid records (plastomes ~120-220 kb) are checked against --length2_cpgenes instead
    plastid = pd.Series(
        [
            _record_group(str(o).lower() if pd.notna(o) else "", str(d) if pd.notna(d) else "") == "plastid"
            for o, d in zip(
                df.get("organelle", pd.Series("", index=df.index)), df.get("Definition", pd.Series("", index=df.index))
            )
        ],
        index=df.index,
        dtype=bool,
    )
    cp_lo, cp_hi = parse_interval(config.length2_cpgenes)
    cp_ok = plastid & (df["Length"] >= (cp_lo or 0)) & (df["Length"] <= (cp_hi if cp_hi is not None else float("inf")))
    keep |= cp_ok
    removed.extend(
        _filtered_rows(
            df[~keep],
            "classify round 1",
            lambda r: (
                f"length {int(r['Length'])} bp outside {config.global_length_range} "
                "(global length filter; change with --length_range2_all)"
            ),
        )
    )
    df = df[keep]
    logger.info(f"After global length filter: {len(df)} records")

    df = _filter_dataframe_logged(df, moleculetype, mol_type, length_filter, organelle_filter, removed)
    logger.info(f"After additional filters: {len(df)} records")

    results = []
    conflict_records = []
    for row_dict in df.to_dict("records"):  # itertuples renamed columns with spaces
        res = process_row(row_dict, config)
        row_dict["gene_type"] = res["gene_type"]
        row_dict["match_source"] = res["match_source"]
        row_dict["Original_match"] = res["Original_match"]
        row_dict["Conflict"] = res["conflict"]
        row_dict["Conflict_reason"] = res["Conflict_reason"]
        row_dict["Assignment_reason"] = res["Assignment_reason"]
        if res["conflict"]:
            conflict_records.append(row_dict)
        results.append(row_dict)

    res_df = pd.DataFrame(results) if results else pd.DataFrame(columns=list(df.columns) + DIAG_COLUMNS)
    if "LocusID" in res_df.columns and len(res_df):
        dup = res_df.duplicated(subset=["LocusID"], keep="first")
        removed.extend(_filtered_rows(res_df[dup], "classify round 1", "duplicate LocusID (first occurrence kept)"))
        res_df = res_df[~dup]
    _write_filtered(output_dir, removed, append=False)
    if removed:
        logger.info(f"Filtered out {len(removed)} records -> {FILTERED_FILE} (with reasons)")
    res_df = res_df.drop(columns=["g2t_input_row", "g2t_length_raw"], errors="ignore")

    if not config.add_assignment_reasons and "Assignment_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Assignment_reason"])
    if not config.record_conflict_reasons and "Conflict_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Conflict_reason"])

    res_df = reorder_columns(res_df)

    assigned_df = res_df[res_df["gene_type"] != ""]
    unmatched_df = res_df[res_df["gene_type"] == ""]

    _save_csv(assigned_df, output_dir, prematch_output, "Matched")
    _save_csv(unmatched_df, output_dir, unmatched_output, "Unmatched")

    if conflict_records:
        conflict_df = reorder_columns(
            pd.DataFrame(conflict_records).drop(columns=["g2t_input_row", "g2t_length_raw"], errors="ignore")
        )
        _save_csv(conflict_df, output_dir, conflicted_output, "Conflict")

    multip_df = assigned_df[assigned_df["gene_type"].apply(lambda x: len(x.split(",")) > 1 if x else False)]
    if not multip_df.empty:
        _save_csv(multip_df, output_dir, multip_output, "Multi-match")

    total = len(res_df)
    logger.info(
        f"Round 1 summary: total={total}, matched={len(assigned_df)}, "
        f"conflict={len(conflict_records)}, unmatched={len(unmatched_df)}"
    )
    return res_df


def process_recheck(
    input_file: str,
    output_dir: str,
    config: MatchConfig,
    apply_length_filter: bool = True,
    assigned_out: str = "assigned_genes_types2.csv",
    multip_out: str = "multiple_matched2.csv",
    unmatched_out: str = "unmatched_sequences2.csv",
) -> pd.DataFrame:
    os.makedirs(output_dir, exist_ok=True)

    ext = os.path.splitext(input_file)[1].lower()
    if ext not in (".csv", ".tsv", ".xls", ".xlsx", ".txt"):
        logger.error(f"Unsupported format: {ext}")
        raise UnsupportedFormatError(f"Unsupported format: {ext}")
    df = utils_read_table(input_file)
    df.columns = df.columns.str.strip()

    if df.empty:
        logger.info("Round 2: no unmatched records to recheck")
        _save_csv(df, output_dir, assigned_out, "Round 2 matched")
        _save_csv(df, output_dir, unmatched_out, "Round 2 unmatched")
        return df

    # Use union of mt+nt ranges instead of nuclear-only
    if apply_length_filter:
        df["Length"] = df["Length"].astype(str).str.replace(r"\s*bp\s*", "", regex=True).str.replace(",", "")
        df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
        ranges = [parse_interval(r) for r in (config.length2_mtgenes, config.length2_ntgenes, config.length2_cpgenes)]
        lo = min(r[0] or 0 for r in ranges)
        hi = max(r[1] or float("inf") for r in ranges)
        if lo is not None:
            df = df[df["Length"] >= lo]
        if hi is not None and hi != float("inf"):
            df = df[df["Length"] <= hi]

    results = []
    conflict_records = []
    for row_dict in df.to_dict("records"):  # itertuples renamed columns with spaces
        res = recheck_match_row(row_dict, config)
        row_dict["gene_type"] = res["gene_type"]
        row_dict["match_source"] = res["match_source"]
        row_dict["Original_match"] = res["Original_match"]
        row_dict["Conflict"] = res["conflict"]
        row_dict["Conflict_reason"] = res["Conflict_reason"]
        row_dict["Assignment_reason"] = res["Assignment_reason"]
        if res["conflict"]:
            conflict_records.append(row_dict)
        results.append(row_dict)

    res_df = pd.DataFrame(results) if results else pd.DataFrame(columns=list(df.columns) + DIAG_COLUMNS)
    if "LocusID" in res_df.columns:
        res_df = res_df.drop_duplicates(subset=["LocusID"])

    if not config.add_assignment_reasons and "Assignment_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Assignment_reason"])
    if not config.record_conflict_reasons and "Conflict_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Conflict_reason"])

    res_df = reorder_columns(res_df)

    assigned_df = res_df[res_df["gene_type"] != ""]
    unmatched_df = res_df[res_df["gene_type"] == ""]
    multi_df = assigned_df[assigned_df["gene_type"].apply(lambda x: len(x.split(",")) > 1 if x else False)]

    _save_csv(assigned_df, output_dir, assigned_out, "Round 2 matched")
    _save_csv(unmatched_df, output_dir, unmatched_out, "Round 2 unmatched")
    if not multi_df.empty:
        _save_csv(multi_df, output_dir, multip_out, "Round 2 multi-match")

    if conflict_records:
        conflict_file = os.path.join(output_dir, "conflicted_genes.csv")
        conflict_df = reorder_columns(pd.DataFrame(conflict_records))
        if os.path.exists(conflict_file):
            try:
                existing = pd.read_csv(conflict_file, dtype=str)
                combined = pd.concat([existing, conflict_df]).drop_duplicates()
                combined.to_csv(conflict_file, index=False)
                logger.info(f"Appended {len(conflict_df)} conflict records to {conflict_file}")
            except Exception as e:
                logger.error(f"Failed to merge conflict records: {e}")
        else:
            _save_csv(conflict_df, output_dir, "conflicted_genes.csv", "Round 2 conflict")

    logger.info(
        f"Round 2 summary: total={len(res_df)}, matched={len(assigned_df)}, "
        f"conflict={len(conflict_records)}, unmatched={len(unmatched_df)}"
    )
    return res_df
