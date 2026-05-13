#!/usr/bin/env python3
"""
Species-Gene Organizer (refactored)
===================================
Groups gene assignment results by species_voucher, one row per species,
with gene-type LocusID columns and metadata summaries.

Improvements:
  - Fixed filenames (organized, single .py extension)
  - Shared utility functions imported from g2t.utils
  - Raises exceptions instead of sys.exit() for testability
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Set, Optional, Any
import logging
import time

from g2t.utils import read_table, write_csv, StepResult

logger = logging.getLogger(__name__)


@dataclass
class OrganizeConfig:
    group_column: str = "species_voucher_new"
    first_col_name: str = "species_voucher_new"
    second_col_name: str = "organism"

    gene_order: List[str] = field(default_factory=lambda: [
        "mtgenome", "coi", "16s", "12s", "cob", "cox2", "cox3",
        "18s", "28s", "its1-its2", "18-28s", "ef-1", "h3"
    ])

    metadata_mode: str = "first_nonempty"

    meta_columns: List[str] = field(default_factory=lambda: [
        "Conflict", "match_source", "Original_match", "Conflict_reason",
        "Assignment_reason", "Class", "Order", "Family", "Genus",
        "PCR_primers", "geo_loc_name", "lat_lon", "collected_by",
        "identified_by", "altitude",
        "country", "collection_date", "isolate", "specimen_voucher",
        "host", "habitat", "isolation_source", "culture_collection",
        "clone", "strain", "note",
    ])

    extra_columns: List[str] = field(default_factory=list)

    mtgenome_prefer_nc: bool = False
    organize_18_28s: bool = True
    organize_mtgenome: bool = True

    gene_includes: Dict[str, List[str]] = field(default_factory=lambda: {
        "18-28s": ["18s", "28s", "its1-its2"],
        "mtgenome": ["coi", "16s", "12s", "cob", "cox2", "cox3"],
    })

    tail_columns: List[str] = field(default_factory=lambda: [
        "TaxonID", "geo_loc_name", "lat_lon",
        "Ref1Authors", "Ref1Title", "Ref1Journal"
    ])

    def get_mito_genes(self) -> List[str]:
        known_mito = {"coi", "16s", "12s", "cob", "cox2", "cox3",
                      "nd1", "nd2", "nd3", "nd4", "nd4l", "nd5", "nd6",
                      "atp6", "atp8", "cytb"}
        return [g for g in self.gene_order if g in known_mito]


def check_columns(df: pd.DataFrame, mandatory: List[str],
                   optional: List[str] = None) -> bool:
    missing = [c for c in mandatory if c not in df.columns]
    if missing:
        logger.error(f"Missing required columns: {', '.join(missing)}")
        return False
    if optional:
        for c in optional:
            if c not in df.columns:
                df[c] = ""
    return True


def get_first_nonempty(series: pd.Series) -> str:
    for v in series:
        if pd.notna(v) and str(v).strip():
            return str(v).strip()
    return ""


def normalize_gene_name(gene: str) -> str:
    return gene.rstrip('#?')


def collect_locusids_per_gene(group: pd.DataFrame, config: OrganizeConfig) -> Dict[str, str]:
    gene_order = config.gene_order
    gene_ids: Dict[str, List[str]] = {g: [] for g in gene_order}

    for row in group.itertuples(index=False):
        if pd.isna(row.gene_type) or pd.isna(row.LocusID):
            continue
        locusid = str(row.LocusID).strip()
        if not locusid:
            continue
        raw_genes = [g.strip() for g in str(row.gene_type).split(',')]
        genes = [normalize_gene_name(g) for g in raw_genes]
        for gene in genes:
            if gene in gene_ids:
                gene_ids[gene].append(locusid)

    for gene in gene_order:
        gene_ids[gene] = list(dict.fromkeys(gene_ids[gene]))

    if config.mtgenome_prefer_nc and len(gene_ids.get("mtgenome", [])) > 1:
        nc_ids = [lid for lid in gene_ids["mtgenome"] if lid.upper().startswith("NC_")]
        if nc_ids:
            gene_ids["mtgenome"] = nc_ids

    if config.organize_18_28s:
        _propagate_locusids(gene_ids, "18-28s", config.gene_includes.get("18-28s", []))
    if config.organize_mtgenome:
        mito_children = config.get_mito_genes()
        _propagate_locusids(gene_ids, "mtgenome", mito_children)

    return {gene: ";".join(ids) for gene, ids in gene_ids.items()}


def _propagate_locusids(gene_ids: Dict[str, List[str]], parent: str,
                        children: List[str]) -> None:
    parent_ids = gene_ids.get(parent, [])
    if not parent_ids:
        return
    parent_set = set(parent_ids)
    for child in children:
        if child in gene_ids:
            child_set = set(gene_ids[child])
            combined = child_set | parent_set
            combined.discard("")
            # Preserve insertion order: child's original order first, then parent extras
            original = list(dict.fromkeys(gene_ids[child]))
            extras = [x for x in parent_ids if x not in child_set]
            gene_ids[child] = list(dict.fromkeys(original + extras))


def _format_gene_values(gene_vals: Dict[str, Set[str]]) -> str:
    active = {g: vs for g, vs in gene_vals.items() if vs}
    if not active:
        return ""

    all_vals = set()
    for vs in active.values():
        all_vals.update(vs)
    if not all_vals:
        return ""

    if len(all_vals) == 1:
        return next(iter(all_vals))

    val_sets = list(active.values())
    if all(vs == val_sets[0] for vs in val_sets):
        return ";".join(sorted(val_sets[0]))

    all_active_genes = set(active.keys())
    parts = []
    for val in sorted(all_vals):
        genes_with = sorted(g for g, vs in active.items() if val in vs)
        if set(genes_with) == all_active_genes:
            parts.append(val)
        else:
            parts.append(f"{','.join(genes_with)}:{val}")

    return ";".join(parts)


def metadata_first_nonempty(group: pd.DataFrame, columns: List[str]) -> Dict[str, str]:
    return {col: get_first_nonempty(group[col]) if col in group.columns else ""
            for col in columns}


def metadata_all(group: pd.DataFrame, columns: List[str]) -> Dict[str, str]:
    result = {}
    for col in columns:
        if col not in group.columns:
            result[col] = ""
            continue
        vals = set()
        for v in group[col]:
            if pd.notna(v) and str(v).strip():
                for part in str(v).split(";"):
                    part = part.strip()
                    if part:
                        vals.add(part)
        result[col] = ";".join(sorted(vals))
    return result


def metadata_per_gene(group: pd.DataFrame, gene_order: List[str],
                      columns: List[str]) -> Dict[str, str]:
    result = {}
    for col in columns:
        if col not in group.columns:
            result[col] = ""
            continue
        gene_vals: Dict[str, Set[str]] = {g: set() for g in gene_order}
        for row in group.itertuples(index=False):
            if pd.isna(row.gene_type):
                continue
            val = getattr(row, col, None)
            if pd.isna(val) or not str(val).strip():
                continue
            val_str = str(val).strip()
            genes = [normalize_gene_name(g.strip()) for g in str(row.gene_type).split(',')]
            for gene in genes:
                if gene in gene_vals:
                    gene_vals[gene].add(val_str)
        result[col] = _format_gene_values(gene_vals)
    return result


def metadata_advanced(group: pd.DataFrame, gene_locusids: Dict[str, str],
                      columns: List[str]) -> Dict[str, str]:
    result = {}
    col_locusid_val: Dict[str, Dict[str, str]] = {col: {} for col in columns}
    for row in group.itertuples(index=False):
        locusid = str(row.LocusID).strip() if pd.notna(row.LocusID) else ""
        if not locusid:
            continue
        for col in columns:
            val = getattr(row, col, None)
            if pd.isna(val) or not str(val).strip():
                continue
            val_str = str(val).strip()
            if locusid not in col_locusid_val[col]:
                col_locusid_val[col][locusid] = val_str
            elif val_str not in col_locusid_val[col][locusid]:
                col_locusid_val[col][locusid] += f";{val_str}"

    for col in columns:
        this_col_data = col_locusid_val.get(col, {})
        if not this_col_data:
            result[col] = ""
            continue

        gene_value_map: Dict[str, Set[str]] = {}
        for gene, locusids_str in gene_locusids.items():
            if not locusids_str:
                continue
            vals_for_gene: Set[str] = set()
            for lid in locusids_str.split(";"):
                lid = lid.strip()
                if lid in this_col_data:
                    vals_for_gene.update(
                        v.strip() for v in this_col_data[lid].split(";") if v.strip()
                    )
            if vals_for_gene:
                gene_value_map[gene] = vals_for_gene

        mapped_lids = set()
        for locusids_str in gene_locusids.values():
            if locusids_str:
                mapped_lids.update(lid.strip() for lid in locusids_str.split(";"))
        orphan_vals: Set[str] = set()
        for lid, val_str in this_col_data.items():
            if lid not in mapped_lids:
                orphan_vals.update(v.strip() for v in val_str.split(";") if v.strip())
        if orphan_vals:
            gene_value_map["_unlinked"] = orphan_vals

        formatted = _format_gene_values(gene_value_map)
        result[col] = formatted.replace("_unlinked:", "")

    return result


def process_group(group: pd.DataFrame, config: OrganizeConfig) -> Dict[str, Any]:
    out: Dict[str, Any] = {}

    out[config.first_col_name] = group[config.group_column].iloc[0]

    if config.second_col_name in group.columns:
        out[config.second_col_name] = get_first_nonempty(group[config.second_col_name])
    else:
        out[config.second_col_name] = ""

    gene_locusids = collect_locusids_per_gene(group, config)
    for gene in config.gene_order:
        out[gene] = gene_locusids.get(gene, "")

    columns_to_process = list(config.meta_columns)
    if config.metadata_mode == "all":
        exclude = {config.group_column, config.second_col_name,
                   "gene_type", "LocusID"} | set(config.gene_order)
        all_cols = [c for c in group.columns if c not in exclude]
        columns_to_process = all_cols

    existing_cols = [c for c in columns_to_process if c in group.columns]

    if config.metadata_mode == "all":
        meta_out = metadata_all(group, existing_cols)
    elif config.metadata_mode == "advanced":
        meta_out = metadata_advanced(group, gene_locusids, existing_cols)
    elif config.metadata_mode == "per_gene":
        meta_out = metadata_per_gene(group, config.gene_order, existing_cols)
    else:
        meta_out = metadata_first_nonempty(group, existing_cols)

    for col in columns_to_process:
        out[col] = meta_out.get(col, "")

    for col in config.tail_columns:
        if col not in out:
            if col in group.columns:
                out[col] = get_first_nonempty(group[col])
            else:
                out[col] = ""

    for col in config.extra_columns:
        if col not in out:
            if col in group.columns:
                out[col] = get_first_nonempty(group[col])
            else:
                out[col] = ""

    return out


def build_header(config: OrganizeConfig, sample_row: Dict[str, Any]) -> List[str]:
    header = [config.first_col_name, config.second_col_name]
    header.extend(config.gene_order)

    if config.metadata_mode == "all":
        exclude = set(header)
        meta_cols = [k for k in sample_row if k not in exclude
                     and k not in config.tail_columns
                     and k not in config.extra_columns]
        header.extend(meta_cols)
    else:
        header.extend(config.meta_columns)

    header.extend(config.tail_columns)
    header.extend(config.extra_columns)
    return list(dict.fromkeys(header))


def organize(input_file: str, output_file: str, config: OrganizeConfig = None) -> StepResult:
    """Run Step 4: organize genes by species."""
    if config is None:
        config = OrganizeConfig()

    start_time = time.time()

    logger.info(f"Input: {input_file}")
    logger.info(f"Group column: {config.group_column}")
    logger.info(f"Metadata mode: {config.metadata_mode}")

    df = read_table(input_file)
    logger.info(f"Loaded {len(df)} records, {len(df.columns)} columns")

    mandatory = [config.group_column, config.second_col_name, "gene_type", "LocusID"]
    optional = config.meta_columns + config.tail_columns + config.extra_columns
    if not check_columns(df, mandatory, optional):
        return StepResult(success=False, output_file=output_file, elapsed=time.time() - start_time)

    nan_count = df[config.group_column].isna().sum()
    if nan_count > 0:
        logger.warning(f"Group column '{config.group_column}' has {nan_count} NaN rows -> 'UNKNOWN'")
        df[config.group_column] = df[config.group_column].fillna("UNKNOWN")

    groups = df.groupby(config.group_column, sort=True)
    logger.info(f"Total groups: {groups.ngroups}")

    rows = []
    for _, group_df in groups:
        rows.append(process_group(group_df, config))

    if not rows:
        logger.error("No output rows generated")
        return StepResult(success=False, output_file=output_file, elapsed=time.time() - start_time)

    header = build_header(config, rows[0])
    out_df = pd.DataFrame(rows, columns=header)

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    out_df.to_csv(output_file, index=False)
    logger.info(f"Output: {output_file} ({len(out_df)} rows, {len(header)} columns)")

    gene_stats = []
    for gene in config.gene_order:
        if gene in out_df.columns:
            non_empty = out_df[gene].astype(str).apply(lambda x: x.strip() != "").sum()
            gene_stats.append(f"{gene}={non_empty}")
    logger.info(f"Gene coverage: {', '.join(gene_stats)}")

    elapsed = time.time() - start_time
    return StepResult(success=True, output_file=output_file, rows=len(out_df), elapsed=elapsed)


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Species-Gene Organizer: group gene assignments by species",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--input", required=True, help="Input file (CSV/TSV/Excel)")
    parser.add_argument("-o", "--output", default="organized_species_voucher.csv",
                        help="Output filename")

    parser.add_argument("--group_column", "--which_column_as_species_voucher",
                        default="species_voucher_new",
                        help="Column to group rows by")
    parser.add_argument("--first_column_name", "--first_row_name",
                        default="species_voucher_new",
                        help="Name for the first column in output")
    parser.add_argument("--second_column_name", "--second_row_name",
                        default="organism",
                        help="Name for the second column in output")

    parser.add_argument("--gene_order", "--rank_genes_order_in_rows",
                        nargs="+",
                        default=["mtgenome", "coi", "16s", "12s", "cob", "cox2", "cox3",
                                 "18s", "28s", "its1-its2", "18-28s", "ef-1", "h3"],
                        help="Gene types to include as columns, in order")

    parser.add_argument("--metadata_mode", "--write_metadata_mode",
                        choices=["first_nonempty", "casual", "advanced", "per_gene", "all"],
                        default="first_nonempty",
                        help="How to aggregate metadata across rows")

    # Hidden backward-compatible aliases (legacy yes/no flags)
    parser.add_argument("--write_metadata_first_nonempty", choices=["yes", "no"], default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--write_metadata_advanced", choices=["yes", "no"], default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--if_organize_all", choices=["yes", "no"], default=None,
                        help=argparse.SUPPRESS)

    parser.add_argument("--metadata_columns", "--mustbe_meta_data_list",
                        nargs="+",
                        default=["Conflict", "match_source", "Original_match",
                                 "Conflict_reason", "Assignment_reason",
                                 "Class", "Order", "Family", "Genus",
                                 "PCR_primers", "geo_loc_name", "lat_lon",
                                 "collected_by", "identified_by", "altitude"],
                        help="Metadata columns to include in output")
    parser.add_argument("--extra_columns", "--extra_information_write_in",
                        nargs="*", default=[],
                        help="Additional columns to append to output")

    parser.add_argument("--mtgenome_prefer_nc", "--mtgenome_dupli_write_nc",
                        choices=["yes", "no"], default="no",
                        help="Prefer NC_ accessions when mtgenome has duplicates")
    parser.add_argument("--organize_18_28s_separately", "--if_organize_18_28s_into_seprate",
                        choices=["yes", "no"], default="yes",
                        help="Propagate 18-28s LocusIDs to sub-gene columns")
    parser.add_argument("--organize_mtgenome_separately", "--if_organize_mtgenomes_into_sep",
                        choices=["yes", "no"], default="yes",
                        help="Propagate mtgenome LocusIDs to sub-gene columns")
    parser.add_argument("--mtgenome_includes", nargs="*", default=None,
                        help="Override which genes mtgenome includes")

    args = parser.parse_args(argv)

    _DEPRECATED_FLAGS = {
        "--which_column_as_species_voucher": "--group_column",
        "--first_row_name": "--first_column_name",
        "--second_row_name": "--second_column_name",
        "--rank_genes_order_in_rows": "--gene_order",
        "--write_metadata_mode": "--metadata_mode",
        "--mustbe_meta_data_list": "--metadata_columns",
        "--extra_information_write_in": "--extra_columns",
        "--mtgenome_dupli_write_nc": "--mtgenome_prefer_nc",
        "--if_organize_18_28s_into_seprate": "--organize_18_28s_separately",
        "--if_organize_mtgenomes_into_sep": "--organize_mtgenome_separately",
    }
    check_argv = argv if argv else sys.argv
    for old, new in _DEPRECATED_FLAGS.items():
        if old in check_argv:
            warnings.warn(f"{old} is deprecated, use {new} instead", DeprecationWarning, stacklevel=2)

    # Backward-compatible mode mapping
    if args.if_organize_all == "yes":
        metadata_mode = "all"
    elif args.write_metadata_advanced == "yes":
        metadata_mode = "advanced"
    elif args.write_metadata_first_nonempty == "no":
        metadata_mode = "per_gene"
    elif args.metadata_mode:
        metadata_mode = "first_nonempty" if args.metadata_mode == "casual" else args.metadata_mode
    else:
        metadata_mode = "first_nonempty"

    gene_includes = {
        "18-28s": ["18s", "28s", "its1-its2"],
        "mtgenome": ["coi", "16s", "12s", "cob", "cox2", "cox3"],
    }
    if args.mtgenome_includes:
        gene_includes["mtgenome"] = [g.strip().lower() for g in args.mtgenome_includes]

    config = OrganizeConfig(
        group_column=args.group_column,
        first_col_name=args.first_column_name,
        second_col_name=args.second_column_name,
        gene_order=args.gene_order,
        metadata_mode=metadata_mode,
        meta_columns=args.metadata_columns,
        extra_columns=args.extra_columns or [],
        mtgenome_prefer_nc=(args.mtgenome_prefer_nc.lower() == "yes"),
        organize_18_28s=(args.organize_18_28s_separately.lower() == "yes"),
        organize_mtgenome=(args.organize_mtgenome_separately.lower() == "yes"),
        gene_includes=gene_includes,
    )

    result = organize(args.input, args.output, config)

    if not result.success:
        sys.exit(1)


if __name__ == "__main__":
    main()
