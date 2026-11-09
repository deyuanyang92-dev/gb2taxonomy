"""Gene type classifier (g2t step 2): GenBank metadata table -> one gene type per record.

This module holds the step (``classify()``), its command line (``g2t-classify``) and logging; the parts
are in ``g2t/_classify/`` and their names are re-exported here:

- ``genedict``  packaged synonyms (data/gene_dict.yaml), user dictionaries, keyword patterns
- ``config``    MatchConfig (gene sets, dictionary, length limits of one run)
- ``text``      phrase matching, DEFINITION cleaning, lengths, negative context
- ``tables``    record filters, output file names, column order
- ``match``     round 1: DEFINITION -> gene types per compartment
- ``conflict``  round 1.5: check against title and length, reassign
- ``recheck``   round 2: unmatched records against all compartments
- ``report``    per-gene files, unmatched report, record_status.csv
- ``flows``     run the rounds over a table
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings

import pandas as pd

from g2t._classify import tables as _tables
from g2t._classify.config import (  # noqa: F401 - re-exported
    MatchConfig,
)
from g2t._classify.conflict import (  # noqa: F401 - re-exported
    conflict_check_and_reassign,
)
from g2t._classify.flows import (  # noqa: F401 - re-exported
    process_prematch,
    process_recheck,
    process_row,
)
from g2t._classify.genedict import (  # noqa: F401 - re-exported
    default_gene_dict,
    init_gene_groups,
    load_path_dict,
    read_gene_dict_file,
)
from g2t._classify.match import (  # noqa: F401 - re-exported
    first_round_match,
    match_mito_genes,
    match_nuclear_genes,
    match_plastid_genes,
)
from g2t._classify.recheck import (  # noqa: F401 - re-exported
    recheck_match_row,
)
from g2t._classify.report import (  # noqa: F401 - re-exported
    build_record_status,
    categorize_unmatched_definition,
    extract_genes,
    generate_unmatched_report,
)
from g2t._classify.tables import (  # noqa: F401 - re-exported
    CLASSIFY_OUTPUTS,
    DIAG_COLUMNS,
    FILTERED_FILE,
    STATUS_FILE,
    filter_dataframe,
    reorder_columns,
)
from g2t._classify.text import (  # noqa: F401 - re-exported
    NEGATIVE_CONTEXT,
    any_match,
    clean_definition,
    clean_length,
    smart_match,
)
from g2t.genes import MITO_TYPES, NUCLEAR_TYPES, PLASTID_TYPES
from g2t.utils import StepResult


# =========================
# Logging configuration
# =========================
class CustomFormatter(logging.Formatter):
    def format(self, record):
        if record.levelno == logging.DEBUG:
            record.levelname = "CHECKING"
        return super().format(record)


def setup_logging(log_file: str = "classify.log") -> logging.Handler:
    """Write the full (DEBUG) log of this run to ``log_file``. The console is left to the caller
    (``g2t-classify`` configures it once; adding a second console handler printed every line twice)."""
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(CustomFormatter("%(asctime)s - %(levelname)s - %(message)s"))
    pkg_logger = logging.getLogger("g2t")
    pkg_logger.setLevel(logging.DEBUG)
    pkg_logger.addHandler(handler)
    return handler


# The documented interface. The other names above stay importable from here for existing code; new code
# that needs internals imports them from g2t._classify.<module>.
__all__ = [
    "MatchConfig",
    "classify",
    "default_gene_dict",
    "load_path_dict",
    "main",
    "read_gene_dict_file",
    "setup_logging",
]

logger = logging.getLogger("g2t.classify")  # also under `python -m g2t.classify` (__name__ is "__main__" there)


# =========================
# Main
# =========================
def _run_cli(args: argparse.Namespace) -> None:
    os.makedirs(args.output, exist_ok=True)
    pkg_logger = logging.getLogger("g2t")
    old_level = pkg_logger.level
    handler = setup_logging(os.path.join(args.output, "classify.log"))
    try:
        _run_cli_steps(args)
    finally:
        pkg_logger.removeHandler(handler)
        pkg_logger.setLevel(old_level)
        handler.close()


def _run_cli_steps(args: argparse.Namespace) -> None:
    logger.info("=" * 60)
    logger.info("Gene Type Classifier (refactored)")
    logger.info("=" * 60)

    output_dir = args.output

    mtgenes = [s.strip() for s in args.mtgenes_list.split(",") if s.strip()]
    ntgenes = [s.strip() for s in args.ntgenes_list.split(",") if s.strip()]
    cpgenes = [s.strip() for s in args.cpgenes_list.split(",") if s.strip()]
    gene_dict = load_path_dict(args.path_dict) if args.path_dict else default_gene_dict()
    mito_genes, nuclear_genes, plastid_genes = init_gene_groups(mtgenes, ntgenes, cpgenes, gene_dict)

    gene_ranges = {
        "mtgenome": args.length_range2_mtgenome,
        "plastome": args.length_range2_plastome,
        "18-28s": args.length_range2_18_28s,
        "coi": args.length_range2_coi,
        "its1-its2": args.length_range2_its,
        "cob": args.length_range2_cob,
        "12s": args.length_range2_12s,
        "16s": args.length_range2_16s,
        "18s": args.length_range2_18s,
        "28s": args.length_range2_28s,
        "h3": args.length_range2_h3,
    }
    if args.gene_length_range:
        for item in args.gene_length_range:
            if "=" in item:
                gene, rng = item.split("=", 1)
                gene_ranges[gene.strip().lower()] = rng.strip()

    config = MatchConfig(
        if_remove_matched=args.if_remove_matched_sequences,
        length2_mtgenes=args.length2_mtgenes,
        length2_ntgenes=args.length2_ntgenes,
        length2_cpgenes=args.length2_cpgenes,
        global_length_range=args.length_range2_all,
        gene_length_ranges=gene_ranges,
        add_assignment_reasons=args.if_add_assigned_gene_types_reasons.lower() == "yes",
        record_conflict_reasons=args.if_recorded_conflict_reasons.lower() == "yes",
        mito_genes=mito_genes,
        nuclear_genes=nuclear_genes,
        plastid_genes=plastid_genes,
        gene_dict=gene_dict,
    )

    gene_types_arg = args.which_gene_types_extract
    result = classify(
        args.input,
        output_dir,
        config,
        prematch_output=args.prematch_output_file_name,
        unmatched_output=args.unmatched_sequences,
        if_recheck=args.if_recheck_unmatched.lower() == "yes",
        apply_length_filter_recheck=args.if_length_range_work2_recheck.lower() == "yes",
        extract_gene_types=[s.strip().lower() for s in gene_types_arg.split(",")] if gene_types_arg else None,
        multiple_output=args.multip_matched_genestype,
        conflicted_output=args.conflicted_sequences,
        moleculetype=args.moleculetype,
        mol_type=args.mol_type,
        length=args.length,
        organelle=args.organelle,
    )
    if not result.success:
        sys.exit(1)
    logger.info("=" * 60)
    logger.info("Processing complete")
    logger.info("=" * 60)


def classify(
    input_file: str,
    output_dir: str,
    config: MatchConfig | None = None,
    prematch_output: str = "assigned_genes_type.csv",
    unmatched_output: str = "unmatched_sequences.csv",
    if_recheck: bool = True,
    apply_length_filter_recheck: bool = True,
    extract_gene_types: list | None = None,
    *,
    multiple_output: str = "multiple_matched.csv",
    conflicted_output: str = "conflicted_genes.csv",
    moleculetype: str = "all",
    mol_type: str = "all",
    length: str = "all",
    organelle: str = "all",
) -> StepResult:
    """Classify gene types from metadata CSV. Returns StepResult. ``g2t-classify`` calls this too.

    The keyword-only arguments are the file names and record filters of ``g2t-classify``.
    """
    import time

    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    if config is None:
        config = MatchConfig()
    _tables._clean_outputs(output_dir, [prematch_output, unmatched_output, multiple_output, conflicted_output])

    try:
        logger.info("-" * 40)
        logger.info("Round 1: Gene type matching")
        logger.info("-" * 40)
        prematch_df = process_prematch(
            input_file,
            output_dir,
            config,
            prematch_output,
            multiple_output,
            unmatched_output,
            conflicted_output,
            moleculetype,
            mol_type,
            length,
            organelle,
        )
        if extract_gene_types:
            extract_genes(prematch_df, extract_gene_types, output_dir)

        unmatched_file = os.path.join(output_dir, unmatched_output)
        if if_recheck:
            if os.path.exists(unmatched_file) and os.path.getsize(unmatched_file) > 50:
                logger.info("-" * 40)
                logger.info("Round 2: Unmatched record recheck")
                logger.info("-" * 40)
                process_recheck(
                    unmatched_file,
                    output_dir,
                    config,
                    apply_length_filter=apply_length_filter_recheck,
                )
            else:
                logger.warning(f"Unmatched file not found or empty: {unmatched_file}, skipping round 2")
        all_file = _tables._write_all(output_dir, prematch_output)
        logger.info(f"Merged results: {all_file}")

        try:
            build_record_status(input_file, output_dir, if_recheck, config, unmatched_output, prematch_output)
        except Exception as e:  # noqa: BLE001 - accounting must not break classification
            logger.warning(f"Failed to write {STATUS_FILE}: {e}")

        if os.path.exists(unmatched_file):
            try:
                unmatched_df = pd.read_csv(unmatched_file, dtype=str)
                if not unmatched_df.empty:
                    generate_unmatched_report(unmatched_df, output_dir)
            except (OSError, pd.errors.EmptyDataError) as e:
                logger.warning(f"Failed to generate unmatched report: {e}")

        rows = len(pd.read_csv(all_file, dtype=str))
        return StepResult(
            success=True,
            output_file=all_file,
            rows=rows,
            elapsed=time.time() - start_time,
        )
    except Exception as e:
        logger.exception(f"Classify failed: {e}")
        return StepResult(
            success=False, output_file="", elapsed=time.time() - start_time, error=f"{type(e).__name__}: {e}"
        )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Gene Type Classifier: GenBank metadata -> gene_type labels",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-i", "--input", required=True, help="Input file (CSV/TSV/Excel)")
    parser.add_argument("-o", "--output", default="assign_genes_type", help="Output directory")

    parser.add_argument(
        "--prematch_output_file_name", default="assigned_genes_type.csv", help="Round-1 assignments (file name in -o)"
    )
    parser.add_argument(
        "--multip_matched_genestype", default="multiple_matched.csv", help="Records matching several gene types"
    )
    parser.add_argument(
        "--unmatched_sequences", default="unmatched_sequences.csv", help="Records matching no gene type"
    )
    parser.add_argument(
        "--conflicted_sequences", default="conflicted_genes.csv", help="Records whose title contradicts the gene type"
    )

    parser.add_argument(
        "--if_remove_matched_sequences",
        "-ifrms",
        default="no",
        help="yes: stop at the first matching gene type (one gene per record)",
    )

    parser.add_argument("--length_range2_mtgenome", default="3000:500000", help="Length window (bp) of mitogenomes")
    parser.add_argument("--length_range2_plastome", default="20000:300000", help="Length window (bp) of plastomes")
    parser.add_argument(
        "--length_range2_18_28s", default="5000:15000", help="Length window (bp) of whole 18S-ITS-28S records"
    )
    parser.add_argument("--length_range2_coi", default="none", help="Length window (bp) of coi; none = no limit")
    parser.add_argument("--length_range2_its", default="none", help="Length window (bp) of its; none = no limit")
    parser.add_argument("--length_range2_cob", default="none", help="Length window (bp) of cob; none = no limit")
    parser.add_argument("--length_range2_12s", default="none", help="Length window (bp) of 12s; none = no limit")
    parser.add_argument("--length_range2_16s", default="none", help="Length window (bp) of 16s; none = no limit")
    parser.add_argument("--length_range2_18s", default="none", help="Length window (bp) of 18s; none = no limit")
    parser.add_argument("--length_range2_28s", default="none", help="Length window (bp) of 28s; none = no limit")
    parser.add_argument("--length_range2_h3", default="none", help="Length window (bp) of h3; none = no limit")

    parser.add_argument("--gene_length_range", nargs="*", default=None, help="Per-gene length ranges: gene=lower:upper")

    parser.add_argument(
        "--length2_mtgenes", default="150:50000", help="Length (bp) a record needs to be typed as a mitochondrial gene"
    )
    parser.add_argument(
        "--length2_ntgenes", default="150:50000", help="Length (bp) a record needs to be typed as a nuclear gene"
    )
    parser.add_argument(
        "--length2_cpgenes", default="150:300000", help="Length (bp) a record needs to be typed as a plastid gene"
    )
    parser.add_argument(
        "--length_range2_all",
        default="150:50000",
        help="Global length filter of round 1 (plastid records use --length2_cpgenes); others are listed as filtered",
    )
    parser.add_argument("--length", default="all", help="Extra length filter, e.g. 300:2000 (all = none)")
    parser.add_argument("--organelle", default="all", help="Keep only this /organelle value, e.g. mitochondrion")
    parser.add_argument("--moleculetype", default="all", help="Keep only this LOCUS molecule type, e.g. DNA")
    parser.add_argument("--mol_type", default="all", help="Keep only these /mol_type values (comma-separated)")

    parser.add_argument(
        "--mtgenes_list",
        default=",".join(MITO_TYPES),
        help="Mitochondrial gene types searched",
    )
    parser.add_argument("--ntgenes_list", default=",".join(NUCLEAR_TYPES), help="Nuclear gene types searched")
    parser.add_argument(
        "--cpgenes_list",
        default=",".join(PLASTID_TYPES),
        help="Plastid gene types searched",
    )
    parser.add_argument("--path_dict", default=None, help="External gene dictionary file")

    # "--which_gene_types_extact" (misspelt) is still accepted
    parser.add_argument(
        "--which_gene_types_extract",
        "--which_gene_types_extact",
        "-wgte",
        default=None,
        help="Also write one file per gene type for these types (comma-separated)",
    )

    parser.add_argument(
        "--if_recheck_unmatched",
        default="yes",
        help="yes: round 2 re-checks unmatched records against all compartments (circular = genome)",
    )
    parser.add_argument(
        "--if_length_range_work2_recheck", default="yes", help="yes: round 2 applies the length windows too"
    )
    parser.add_argument(
        "--if_add_assigned_gene_types_reasons", default="yes", help="yes: write the Assignment_reason column"
    )
    parser.add_argument("--if_recorded_conflict_reasons", default="yes", help="yes: write the Conflict_reason column")

    args = parser.parse_args(argv)

    check_argv = argv if argv else sys.argv
    if "--which_gene_types_extact" in check_argv:
        warnings.warn(
            "--which_gene_types_extact is deprecated (typo), use --which_gene_types_extract",
            DeprecationWarning,
            stacklevel=2,
        )

    _run_cli(args)


if __name__ == "__main__":
    main()
