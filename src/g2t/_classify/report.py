"""Per-gene files, the unmatched report and record_status.csv."""

from __future__ import annotations

import logging
import os
import re

import pandas as pd

from g2t._classify.config import MatchConfig
from g2t._classify.tables import FILTERED_FILE, STATUS_FILE, reorder_columns
from g2t.utils import read_table as utils_read_table

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


def extract_genes(df: pd.DataFrame, gene_types: list[str] | None, output_dir: str) -> None:
    if not gene_types:
        all_genes: set[str] = set()
        for gt in df["gene_type"].dropna():
            all_genes.update(s.strip() for s in str(gt).split(",") if s.strip())
        gene_types = sorted(all_genes)
        logger.info(f"Auto-detected gene types: {gene_types}")

    for gene in gene_types:
        subset = df[
            df["gene_type"].apply(lambda x, g=gene: g in [s.strip() for s in str(x).split(",")] if x else False)
        ]
        if subset.empty:
            logger.warning(f"No records for gene type '{gene}'")
            continue
        safe_name = re.sub(r"[^\w\-]", "_", gene)
        out_file = os.path.join(output_dir, f"{safe_name}_extracted.csv")
        try:
            reorder_columns(subset).to_csv(out_file, index=False)
            logger.info(f"  {gene} -> {out_file} ({len(subset)} records)")
        except Exception as e:
            logger.error(f"Failed to save {gene} extraction: {e}")


def _save_csv(df: pd.DataFrame, output_dir: str, filename: str, label: str) -> None:
    path = os.path.join(output_dir, filename)
    try:
        df.to_csv(path, index=False)
        logger.info(f"{label}: {path} ({len(df)} records)")
    except Exception as e:
        logger.error(f"Failed to save {label}: {e}")


def categorize_unmatched_definition(def_text: str) -> str:
    """Categorize unmatched sequence by its Definition field."""
    def_lower = str(def_text).lower()

    # Genome/assembly data
    if "whole genome shotgun" in def_lower:
        return "whole_genome_shotgun"
    elif "genome assembly" in def_lower:
        return "genome_assembly"
    elif "tpa_asm" in def_lower or "tpa:" in def_lower:
        return "TPA_assembly"
    elif "chromosome:" in def_lower:
        return "chromosome_data"
    elif "scaffold" in def_lower:
        return "scaffold"
    elif "contig" in def_lower:
        return "contig"

    # RNA types
    elif "mrna sequence" in def_lower or "cdna clone" in def_lower:
        return "mRNA/cDNA"
    elif "trna" in def_lower:
        return "tRNA"
    elif "rrna" in def_lower and "mitochondrial" not in def_lower:
        return "nuclear_rRNA"

    # Other genes not in standard list
    elif "nad" in def_lower or "nadh" in def_lower:
        return "NADH_dehydrogenase"
    elif "atp synthase" in def_lower or "atpase" in def_lower:
        return "ATP_synthase"
    elif "creatine kinase" in def_lower:
        return "creatine_kinase"
    elif "cytochrome b" in def_lower or "cytb" in def_lower:
        return "cytochrome_b_variant"
    else:
        return "other"


def generate_unmatched_report(unmatched_df: pd.DataFrame, output_dir: str) -> str:
    """Generate a summary report of unmatched sequences by category."""
    if unmatched_df.empty:
        return ""

    # Categorize each unmatched sequence
    unmatched_df = unmatched_df.copy()
    unmatched_df["category"] = unmatched_df["Definition"].apply(categorize_unmatched_definition)

    # Generate summary
    summary_lines = []
    summary_lines.append("=" * 60)
    summary_lines.append("UNMATCHED SEQUENCES REPORT")
    summary_lines.append("=" * 60)
    summary_lines.append(f"Total unmatched: {len(unmatched_df)}")
    summary_lines.append("")

    category_counts = unmatched_df["category"].value_counts()
    summary_lines.append("Category breakdown:")
    summary_lines.append("-" * 40)

    for cat, count in category_counts.items():
        pct = count / len(unmatched_df) * 100
        summary_lines.append(f"  {cat}: {count} ({pct:.1f}%)")

    summary_lines.append("")
    summary_lines.append("Sample definitions per category:")
    summary_lines.append("-" * 40)

    for cat in category_counts.index[:5]:  # Top 5 categories
        samples = unmatched_df[unmatched_df["category"] == cat]["Definition"].head(3).tolist()
        summary_lines.append(f"\n[{cat}]")
        for s in samples:
            summary_lines.append(f"  - {str(s)[:100]}...")

    summary_lines.append("")
    summary_lines.append("=" * 60)

    # Save report
    report_path = os.path.join(output_dir, "unmatched_report.txt")
    report_content = "\n".join(summary_lines)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_content)

    # Also save categorized CSV
    categorized_path = os.path.join(output_dir, "unmatched_categorized.csv")
    unmatched_df.to_csv(categorized_path, index=False)

    logger.info(f"Unmatched report saved: {report_path}")
    logger.info(f"Categorized unmatched: {categorized_path}")

    # Print summary to console
    print("\n" + report_content)

    return report_path


def build_record_status(
    input_file: str,
    output_dir: str,
    if_recheck: bool,
    config: MatchConfig,
    unmatched_output: str = "unmatched_sequences.csv",
    prematch_output: str = "assigned_genes_type.csv",
) -> pd.DataFrame:
    """One row per input row: assigned (gene_type) / unmatched / filtered, each with a reason."""
    src = utils_read_table(input_file)
    src.columns = src.columns.str.strip()

    def read(name):
        path = os.path.join(output_dir, name)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            try:
                return pd.read_csv(path, dtype=str).fillna("")
            except pd.errors.EmptyDataError:
                pass
        return pd.DataFrame()

    filt = read(FILTERED_FILE)
    filt_reason = {str(r["input_row"]): r["reason"] for _, r in filt.iterrows()} if len(filt) else {}
    assigned = read("assigned_genes_types_all.csv")
    if assigned.empty:
        assigned = read(prematch_output)
    gene_of = dict(zip(assigned.get("LocusID", []), assigned.get("gene_type", [])))
    reason_of = dict(zip(assigned.get("LocusID", []), assigned.get("Assignment_reason", [""] * len(assigned))))
    recheck_ran = if_recheck and os.path.exists(os.path.join(output_dir, "unmatched_sequences2.csv"))
    un2 = read("unmatched_sequences2.csv") if recheck_ran else pd.DataFrame()
    un2_ids = set(un2.get("LocusID", []))
    rows = []
    for i, r in src.iterrows():
        lid = str(r.get("LocusID", ""))
        base = {
            "input_row": i,
            "ACCESSION": r.get("ACCESSION", ""),
            "LocusID": lid,
            "Organism": r.get("Organism", r.get("organism", "")),
            "Length": r.get("Length", ""),
            "Definition": r.get("Definition", ""),
        }
        if str(i) in filt_reason:
            rows.append({**base, "status": "filtered", "gene_type": "", "reason": filt_reason[str(i)]})
        elif lid in gene_of:
            rows.append({**base, "status": "assigned", "gene_type": gene_of[lid], "reason": reason_of.get(lid, "")})
        else:
            reason = "no gene type matched (DEFINITION keywords, round 1"
            if recheck_ran and lid not in un2_ids:
                reason += (
                    f"); skipped in round-2 recheck: length outside {config.length2_mtgenes} / {config.length2_ntgenes}"
                    f" / {config.length2_cpgenes}"
                )
            else:
                reason += " and 2)" if recheck_ran else ")"
            rows.append({**base, "status": "unmatched", "gene_type": "", "reason": reason})
    out = pd.DataFrame(
        rows,
        columns=[
            "input_row",
            "ACCESSION",
            "LocusID",
            "Organism",
            "Length",
            "Definition",
            "status",
            "gene_type",
            "reason",
        ],
    )
    out.to_csv(os.path.join(output_dir, STATUS_FILE), index=False)
    counts = out["status"].value_counts().to_dict()
    logger.info(
        f"Record status: input={len(out)}, assigned={counts.get('assigned', 0)}, "
        f"unmatched={counts.get('unmatched', 0)}, filtered={counts.get('filtered', 0)} -> {STATUS_FILE}"
    )
    return out
