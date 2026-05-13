#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Species Voucher Builder (refactored)
====================================
Builds unified species_voucher identifiers from gene assignment results.

Improvements:
  - Unicode normalization (sanitize_string)
  - str2bool raises errors for invalid input instead of defaulting to True
  - Fixed output filename typos (vourcher -> voucher, speces -> species)
  - Consolidated frequency group files into single CSV + frequency_group column
  - Shared utility functions imported from g2t.utils
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
import warnings
import pandas as pd
from typing import Dict, List, Optional

from g2t.utils import (
    str2bool, sanitize_string, safe_concat,
    read_table, write_csv, col_key, build_col_lookup, resolve_col,
    StepResult,
)


def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def log(msg: str, quiet: bool, log_file_path: Optional[str] = None) -> None:
    line = f"[{now_ts()}] {msg}"
    if not quiet:
        print(line)
    if log_file_path:
        try:
            with open(log_file_path, "a", encoding="utf-8") as f:
                f.write(line + "\n")
        except Exception:
            pass


def normalize_columns_inplace(df: pd.DataFrame) -> None:
    """Normalize column names in-place: strip whitespace, replace spaces with underscores, deduplicate."""
    new_cols = [col_key(c) for c in df.columns]
    seen: Dict[str, int] = {}
    fixed: List[str] = []
    for c in new_cols:
        if c not in seen:
            seen[c] = 0
            fixed.append(c)
        else:
            seen[c] += 1
            fixed.append(f"{c}__dup{seen[c]}")
    df.columns = fixed


def build_species_vouchers(
    file_path: str,
    output_dir: str,
    if_generate_simple_csv: bool = True,
    if_generate_sorted_group: bool = True,
    combine_order: str = "species_name+pse",
    if_write_haplotype: bool = False,
    columns_order: List[str] = None,
    col_for_empty_fill: str = "ACCESSION",
    primary_voucher_column_name: str = "species_voucher_pse",
    name_updated: str = "species_voucher_new",
    species_name_column: str = "organism",
    simple_csv_columns: str = "organism,gene_type,species_voucher_new,LocusID",
    normalize_column_names: bool = False,
    quiet: bool = False,
    log_file: str = "",
) -> StepResult:
    if columns_order is None:
        columns_order = ["specimen_voucher", "isolate", "culture_collection", "clone", "strain"]
    os.makedirs(output_dir, exist_ok=True)
    log_file_path = log_file

    start_time = time.time()
    log(f"START: input={file_path}", quiet, log_file_path)
    log(f"OUTDIR: {output_dir}", quiet, log_file_path)

    df = read_table(file_path)
    log(f"Loaded table: rows={len(df)} cols={len(df.columns)}", quiet, log_file_path)

    if normalize_column_names:
        normalize_columns_inplace(df)
        log("normalize_column_names=ON: columns spaces -> underscores (and deduped)", quiet, log_file_path)
    else:
        log("normalize_column_names=OFF: keeping original column names", quiet, log_file_path)

    lookup = build_col_lookup(df)

    hap_col = resolve_col(df, "haplotype", lookup)
    species_col = resolve_col(df, species_name_column, lookup)
    fallback_col = resolve_col(df, col_for_empty_fill, lookup)

    resolved_order_cols: List[str] = []
    missing_order_cols: List[str] = []
    for c in columns_order:
        rc = resolve_col(df, c, lookup)
        if rc and rc not in resolved_order_cols:
            resolved_order_cols.append(rc)
        else:
            missing_order_cols.append(c)

    log(f"pse priority columns (requested): {columns_order}", quiet, log_file_path)
    log(f"pse priority columns (found): {resolved_order_cols if resolved_order_cols else '(none)'}", quiet, log_file_path)
    if missing_order_cols:
        log(f"pse priority columns (missing/unresolved): {missing_order_cols}", quiet, log_file_path)

    log(f"species_name_column requested: {species_name_column} -> resolved: {species_col}", quiet, log_file_path)
    log(f"fallback column requested: {col_for_empty_fill} -> resolved: {fallback_col}", quiet, log_file_path)
    log(f"combine_order: {combine_order}", quiet, log_file_path)
    log(f"write haplotype into pse if still empty: {if_write_haplotype}", quiet, log_file_path)

    # 1) only_haplotype.csv
    if hap_col and hap_col in df.columns:
        hap_df = df[df[hap_col].notna()].copy()
        hap_path = os.path.join(output_dir, "only_haplotype.csv")
        write_csv(hap_df, hap_path)
        log(f"only_haplotype.csv written: rows={len(hap_df)} -> {hap_path}", quiet, log_file_path)
    else:
        log("No haplotype column found/resolved; skip only_haplotype.csv", quiet, log_file_path)

    # 2) pse series (vectorized bfill)
    primary_voucher_series = pd.Series(pd.NA, index=df.index)
    fill_from_priority = 0
    fill_from_hap = 0
    fill_from_fallback = 0

    if resolved_order_cols:
        tmp = df[resolved_order_cols].astype("string")
        primary_voucher_series = tmp.bfill(axis=1).iloc[:, 0]
        primary_voucher_series = primary_voucher_series.infer_objects(copy=False)
        fill_from_priority = int(primary_voucher_series.notna().sum())
    else:
        log("WARNING: no priority columns found; pse starts as all-NA", quiet, log_file_path)

    if if_write_haplotype and hap_col and hap_col in df.columns:
        before = primary_voucher_series.notna().sum()
        primary_voucher_series = primary_voucher_series.fillna(df[hap_col].astype("string")).infer_objects(copy=False)
        after = primary_voucher_series.notna().sum()
        fill_from_hap = int(after - before)

    if fallback_col and fallback_col in df.columns:
        before = primary_voucher_series.notna().sum()
        primary_voucher_series = primary_voucher_series.fillna(df[fallback_col].astype("string")).infer_objects(copy=False)
        after = primary_voucher_series.notna().sum()
        fill_from_fallback = int(after - before)

    still_empty = int(primary_voucher_series.isna().sum())
    log(f"pse filled: from_priority={fill_from_priority}, +haplotype={fill_from_hap}, "
        f"+fallback={fill_from_fallback}, still_empty={still_empty}", quiet, log_file_path)

    df[primary_voucher_column_name] = primary_voucher_series

    # 3) updated voucher
    if combine_order == "pse_only":
        df[name_updated] = df[primary_voucher_column_name].apply(sanitize_string)
        log(f"updated voucher = pse_only -> column '{name_updated}' created", quiet, log_file_path)
    else:
        if not species_col or species_col not in df.columns:
            log(f"WARNING: species column unresolved/missing; updated voucher will use pse only", quiet, log_file_path)
            df[name_updated] = df[primary_voucher_column_name].apply(sanitize_string)
        else:
            if combine_order == "species_name+pse":
                df[name_updated] = [
                    safe_concat(species_name, voucher_value)
                    for species_name, voucher_value in zip(df[species_col].tolist(), df[primary_voucher_column_name].tolist())
                ]
            elif combine_order == "pse+species_name":
                df[name_updated] = [
                    safe_concat(voucher_value, species_name)
                    for species_name, voucher_value in zip(df[species_col].tolist(), df[primary_voucher_column_name].tolist())
                ]
            else:
                raise ValueError("combine_order must be: species_name+pse / pse+species_name / pse_only")
            log(f"updated voucher combined -> column '{name_updated}' created", quiet, log_file_path)

    updated_empty = int(df[name_updated].replace("", pd.NA).isna().sum())
    log(f"updated voucher empty rows: {updated_empty}", quiet, log_file_path)

    # 4) full output
    full_output_path = os.path.join(output_dir, "updated_species_voucher.csv")
    write_csv(df, full_output_path)
    log(f"Full output written: {full_output_path}", quiet, log_file_path)

    # 5) simple CSV
    if if_generate_simple_csv:
        want = [c.strip() for c in simple_csv_columns.split(",") if c.strip()]
        resolved_simple: List[str] = []
        missing_simple: List[str] = []
        for c in want:
            rc = resolve_col(df, c, lookup)
            if rc is None and c in df.columns:
                rc = c
            if rc and rc not in resolved_simple:
                resolved_simple.append(rc)
            else:
                if rc is None:
                    missing_simple.append(c)

        if not resolved_simple:
            log("Simple CSV: no requested columns resolved, skipping", quiet, log_file_path)
        else:
            simple_df = df[resolved_simple].copy()
            simple_path = os.path.join(output_dir, "shorten_species_voucher.csv")
            write_csv(simple_df, simple_path)
            log(f"Simple CSV written: cols={len(resolved_simple)} rows={len(simple_df)} -> {simple_path}", quiet, log_file_path)
            if missing_simple:
                log(f"Simple CSV missing columns (unresolved): {missing_simple}", quiet, log_file_path)
    else:
        log("Simple CSV generation OFF", quiet, log_file_path)

    # 6) frequency groups — consolidated into single file with frequency_group column
    if if_generate_sorted_group:
        vc = df[name_updated].replace("", pd.NA).dropna().value_counts()
        group_rows = []
        for freq in range(2, 9):
            keys = vc[vc == freq].index
            if len(keys) == 0:
                continue
            group_df = df[df[name_updated].isin(keys)].copy()
            group_df["frequency_group"] = freq
            group_rows.append(group_df)
            # Still generate individual files for backward compatibility
            outp = os.path.join(output_dir, f"species_voucher_new_group_{freq}.csv")
            write_csv(df[df[name_updated].isin(keys)].copy(), outp)
            log(f"group_{freq}: vouchers={len(keys)} rows={len(group_df)} -> {outp}", quiet, log_file_path)

        if group_rows:
            combined = pd.concat(group_rows, ignore_index=True)
            combined_path = os.path.join(output_dir, "species_voucher_groups.csv")
            write_csv(combined, combined_path)
            log(f"Consolidated groups: {combined_path} ({len(combined)} rows)", quiet, log_file_path)
        else:
            log("No groups (2..8) generated: no vouchers with freq 2..8", quiet, log_file_path)
    else:
        log("Group generation OFF", quiet, log_file_path)

    # 7) count file
    count_series = df[name_updated].replace("", pd.NA).dropna().value_counts()
    freq_distribution = count_series.value_counts().sort_index(ascending=False)
    count_file_path = os.path.join(output_dir, "count_species_voucher.txt")
    with open(count_file_path, "w", encoding="utf-8") as f:
        for freq, count in freq_distribution.items():
            f.write(f"{freq} occurrences: {count}\n")
    log(f"Count file written: {count_file_path}", quiet, log_file_path)

    elapsed = time.time() - start_time
    log(f"DONE. elapsed={elapsed:.2f}s", quiet, log_file_path)

    return StepResult(
        success=True,
        output_file=full_output_path,
        rows=len(df),
        elapsed=elapsed,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Species Voucher Builder: generate primary_voucher and updated voucher identifiers."
    )
    parser.add_argument("-i", "--input", "--file", required=True, help="Input file (xlsx/xls/csv/tsv/txt)")
    parser.add_argument("-o", "--out", required=True, help="Output directory")
    parser.add_argument("--quiet", action="store_true", help="Quiet mode (suppress console output)")
    parser.add_argument("--log_file", default="", help="Log file path")

    parser.add_argument("--normalize_column_names", type=str2bool, default=False,
                        help="Replace column name spaces with underscores (default False)")
    parser.add_argument("--fill_haplotype", "--if_write_haplotype_into_new_species_voucher",
                        type=str2bool, default=True,
                        help="Fill empty primary_voucher with haplotype (default True)")
    parser.add_argument("--columns_to_combine", "--which_columes_want_combine",
                        default="specimen_voucher,isolate,clone,strain",
                        help="Columns used to generate primary_voucher, in priority order")
    parser.add_argument("--fallback_column", "--which_volumn_write2_still_empty_species_voucher",
                        default="ACCESSION",
                        help="Fallback column to fill empty primary_voucher")
    parser.add_argument("--combine_species_name", "--if_species_voucher_combine_species_name",
                        type=str2bool, default=True,
                        help="Combine species name with primary_voucher")
    parser.add_argument("--primary_voucher_column_name", "--name_of_pse_species_voucher",
                        default="species_voucher_pse",
                        help="Column name for the primary voucher identifier")
    parser.add_argument("--updated_voucher_column_name", "--name_of_upated_speces_voucher",
                        default="species_voucher_new",
                        help="Column name for the updated voucher identifier")
    parser.add_argument("--species_name_column", "--which_colum_is_species_name",
                        default="organism",
                        help="Column containing species names")
    parser.add_argument("--simple_csv_columns", "--which_volumns_you_want_for_simple_csv",
                        default="LocusID,Length,Definition,organism,TaxonID,specimen_voucher,species_voucher_new,species_voucher_pse,isolate,strain",
                        help="Columns to include in the simple output CSV")
    parser.add_argument("--generate_simple_csv", "--if_generate_simple_csv",
                        type=str2bool, default=True,
                        help="Generate simple CSV output")
    parser.add_argument("--generate_frequency_groups", "--if_generate_sorted_group",
                        type=str2bool, default=True,
                        help="Generate frequency group files (occurrences 2-8)")
    parser.add_argument("--combine_order", "--set_combine_order_for_species_voucher_new",
                        default="species_name+pse", choices=["species_name+pse", "pse+species_name"],
                        help="How to combine species name with voucher")

    args = parser.parse_args(argv)

    _DEPRECATED_FLAGS = {
        "--which_columes_want_combine": "--columns_to_combine",
        "--which_volumn_write2_still_empty_species_voucher": "--fallback_column",
        "--if_species_voucher_combine_species_name": "--combine_species_name",
        "--name_of_pse_species_voucher": "--primary_voucher_column_name",
        "--name_of_upated_speces_voucher": "--updated_voucher_column_name",
        "--which_colum_is_species_name": "--species_name_column",
        "--which_volumns_you_want_for_simple_csv": "--simple_csv_columns",
        "--if_generate_sorted_group": "--generate_frequency_groups",
        "--set_combine_order_for_species_voucher_new": "--combine_order",
    }
    check_argv = argv if argv else sys.argv
    for old, new in _DEPRECATED_FLAGS.items():
        if old in check_argv:
            warnings.warn(f"{old} is deprecated, use {new} instead", DeprecationWarning, stacklevel=2)

    columns_order = [c.strip() for c in args.columns_to_combine.split(",") if c.strip()]

    if not args.combine_species_name:
        combine_order = "pse_only"
    else:
        combine_order = args.combine_order

    log_file_path = args.log_file.strip()
    if not log_file_path:
        os.makedirs(args.out, exist_ok=True)
        log_file_path = os.path.join(args.out, "update_species_voucher.log")

    build_species_vouchers(
        file_path=args.input,
        output_dir=args.out,
        if_generate_simple_csv=args.generate_simple_csv,
        if_generate_sorted_group=args.generate_frequency_groups,
        combine_order=combine_order,
        if_write_haplotype=args.fill_haplotype,
        columns_order=columns_order,
        col_for_empty_fill=args.fallback_column,
        primary_voucher_column_name=args.primary_voucher_column_name,
        name_updated=args.updated_voucher_column_name,
        species_name_column=args.species_name_column,
        simple_csv_columns=args.simple_csv_columns,
        normalize_column_names=args.normalize_column_names,
        quiet=args.quiet,
        log_file=log_file_path,
    )


if __name__ == "__main__":
    main()
