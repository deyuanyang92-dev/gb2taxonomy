#!/usr/bin/env python3
"""
GenBank metadata extractor (g2t step 1).

Reads GenBank files (.gb, .gbf, .gbk, .gbff, .genbank, also gzip-compressed) and writes one row of metadata per
record: source qualifiers, references, assembly data and taxonomy, with a status for every record.

This module holds the step (``extract()``) and its command line (``g2t-extract``); the parts are in
``g2t/_extract/`` and their names are re-exported here:

- ``files``     finding, opening and checking input files
- ``parse``     one record -> metadata row and assembly row; column order
- ``taxonomy``  lineage from ete3 or Entrez, with a cache (replace ``g2t._extract.taxonomy._taxonomy``,
                not a name here, to swap the service)
- ``status``    per-record / per-file status, summary, reports
- ``process``   reading files (stream or memory) and writing the tables
- ``batch``     one process per file, then merge
"""

from __future__ import annotations

import argparse
import contextlib
import logging
import os
import sys
import time

import pandas as pd

from g2t._extract.batch import (  # noqa: F401 - re-exported
    process_batch_mode,
)
from g2t._extract.files import (  # noqa: F401 - re-exported
    GB_EXTENSIONS,
    copy_failed_files,
    discover_gb_files,
    format_size,
    get_memory_mb,
    is_genbank_name,
    open_genbank,
    psutil,
    validate_and_estimate_records,
)
from g2t._extract.parse import (  # noqa: F401 - re-exported
    ASSEMBLY_CSV_COLUMNS,
    ASSEMBLY_FIELDS,
    BASE_METADATA_COLUMNS,
    canonical_metadata_columns,
    extract_assembly_from_comment_regex,
    extract_assembly_from_record,
    extract_assembly_from_structured_comment,
    extract_metadata,
    get_ordered_fieldnames,
    parse_locus,
    parse_references,
    parse_source_feature,
    prescan_fieldnames,
)
from g2t._extract.process import (  # noqa: F401 - re-exported
    RecordCallback,
    combine_and_save_final_csv,
    count_record_terminators,
    memory_process_single_file,
    process_all_files,
    process_record_unified,
    stream_process_single_file,
)
from g2t._extract.status import (  # noqa: F401 - re-exported
    FileResult,
    FileStatus,
    ProcessingSummary,
    RecordError,
    RecordTrack,
)
from g2t._extract.taxonomy import (  # noqa: F401 - re-exported
    ENTREZ_TAXONOMY_BATCH,
    TAXONOMY_CACHE_NAME,
    TaxonomyService,
    collect_taxon_ids,
    get_taxonomy_lineage,
    initialize_ncbi_taxa,
    prefetch_taxonomy,
    update_ncbi_taxonomy_database,
)
from g2t.utils import StepResult

# ===============================
# Logging
# ===============================
# The documented interface. The other names above stay importable from here for existing code; new code
# that needs internals imports them from g2t._extract.<module>.
__all__ = [
    "GB_EXTENSIONS",
    "TaxonomyService",
    "discover_gb_files",
    "extract",
    "is_genbank_name",
    "main",
    "open_genbank",
]

logger = logging.getLogger("g2t.extract")  # also under `python -m g2t.extract`

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# ===============================
# Public API
# ===============================
def _finish_extract(
    summary: ProcessingSummary,
    output_dir: str,
    start_time: float,
    final_path: str = "",
    rows: int = 0,
    success: bool = False,
    report: bool = True,
    details: bool = False,
    tracking: bool = False,
    error: str = "",
) -> StepResult:
    """Write the reports and build the StepResult.

    ``report``: extraction_report.json (plus the error table when ``details``).
    ``details`` (the command line): copy failed files to failed_extract/ and log the summary.
    ``tracking``: record_tracking.csv.
    """
    summary.total_elapsed_time = time.monotonic() - start_time
    try:
        os.makedirs(output_dir, exist_ok=True)
        if details:
            copy_failed_files(summary.get_failed_files(), output_dir)
            summary.print_summary()
        if report:
            summary.save_report(output_dir)
            if details:
                summary.save_error_csv(output_dir)
        if tracking:
            summary.save_tracking_csv(output_dir)
    except OSError as e:
        logger.error(f"Failed to save extraction report: {e}")
    return StepResult(
        success=success,
        output_file=final_path if success else "",
        rows=rows,
        elapsed=time.monotonic() - start_time,
        error="" if success else error,
    )


def _banner(text: str) -> None:
    logger.info("=" * 75)
    logger.info(text)
    logger.info("=" * 75)


def extract(
    input_files: list,
    output_dir: str,
    stream: bool = True,
    batch: bool = False,
    max_tasks: int = 1,
    final_name: str = "final.csv",
    *,
    metadata_name: str = "metadata.csv",
    assembly_name: str = "assembly.csv",
    metadata_delimiter: str = ",",
    assembly_delimiter: str = "\t",
    final_delimiter: str = ",",
    only_assembly: bool = False,
    only_metadata: bool = False,
    recursive: bool = True,
    log_interval: int = 5000,
    update_taxonomy: bool = False,
    report: bool = True,
    details: bool = False,
    tracking: bool = False,
) -> StepResult:
    """Extract metadata from GenBank files. Returns StepResult. ``g2t-extract`` calls this too.

    Invalid/skipped files are logged with ``logger.warning`` and every file's status,
    record counts and errors are written to ``<output_dir>/extraction_report.json``.
    The step is successful as long as at least one file produced records; it fails
    only when every file failed (or no usable file/final table exists).
    The keyword-only arguments are the options of ``g2t-extract``; their defaults are what the
    pipeline uses (``details``/``tracking`` are on for the command line).
    """
    start_time = time.monotonic()
    os.makedirs(output_dir, exist_ok=True)
    summary = ProcessingSummary()

    def finish(final_path: str = "", rows: int = 0, success: bool = False, error: str = "") -> StepResult:
        return _finish_extract(
            summary,
            output_dir,
            start_time,
            final_path,
            rows,
            success,
            report=report,
            details=details,
            tracking=tracking,
            error=error,
        )

    try:
        if psutil:
            logger.info(f"Available memory: {psutil.virtual_memory().available / (1024**3):.2f} GB")
        _banner("Step 1/4: Discover GenBank files")
        gb_files, total, skipped = discover_gb_files(input_files, recursive)
        summary.total_files_found = total
        summary.non_gb_files_skipped = skipped
        summary.gb_files_found = len(gb_files)
        if not gb_files:
            msg = f"no GenBank files ({', '.join(GB_EXTENSIONS)}, also .gz) found in the inputs"
            logger.error(f"Extract: {msg}")
            return finish(error=msg)

        valid_files = []
        for f in gb_files:
            ok, err, _ = validate_and_estimate_records(f)
            if ok:
                valid_files.append(f)
            else:
                summary.add_result(FileResult(file_path=f, status=FileStatus.INVALID, error_message=err))
                logger.warning(f"Extract: skipping invalid file {os.path.basename(f)} ({err})")
        logger.info(f"Valid files: {len(valid_files)}/{len(gb_files)}")

        if not valid_files:
            msg = "none of the input files is a valid GenBank file"
            logger.error(f"Extract: {msg}")
            return finish(error=msg)

        total_bytes = sum(os.path.getsize(f) for f in valid_files if os.path.exists(f))
        if not stream and total_bytes > 500 * 1024 * 1024:
            logger.warning(f"Large dataset ({format_size(total_bytes)}), recommend --stream")

        # the assembly-only and metadata-only modes have never loaded the taxonomy database
        if not only_assembly and not only_metadata:
            _banner("Step 2/4: Taxonomy module")
            initialize_ncbi_taxa()
            if update_taxonomy:
                update_ncbi_taxonomy_database()

        _banner(f"Step 3/4: Extract metadata ({'stream' if stream else 'memory'} / {'batch' if batch else 'merged'})")
        include_taxonomy = not only_metadata
        if batch:
            processing = process_batch_mode(
                valid_files,
                output_dir,
                max_tasks,
                stream,
                only_assembly,
                only_metadata,
                include_taxonomy,
                metadata_delimiter,
                assembly_delimiter,
                final_delimiter,
                log_interval,
                metadata_output=metadata_name,
                assembly_output=assembly_name,
                final_output=final_name,
            )
        else:
            processing = process_all_files(
                valid_files,
                output_dir,
                metadata_name,
                assembly_name,
                final_name,
                include_taxonomy,
                only_assembly,
                only_metadata,
                metadata_delimiter,
                assembly_delimiter,
                final_delimiter,
                stream,
                log_interval,
            )
        for r in processing.file_results:
            summary.add_result(r)

        _banner("Step 4/4: Report & archive")
        # the merged table needs both halves; the single-half modes produce their own table
        if only_assembly:
            product, sep = os.path.join(output_dir, assembly_name), assembly_delimiter
        elif only_metadata:
            product, sep = os.path.join(output_dir, metadata_name), metadata_delimiter
        else:
            product, sep = os.path.join(output_dir, final_name), final_delimiter
        rows = 0
        if os.path.exists(product):
            with contextlib.suppress(Exception):
                rows = len(pd.read_csv(product, dtype=str, keep_default_na=False, sep=sep))

        n_ok = summary.files_succeeded + summary.files_partial
        problems = [r for r in summary.file_results if r.status != FileStatus.SUCCESS]
        if problems:
            logger.warning(
                f"Extract finished with problems: {summary.files_succeeded} ok, "
                f"{summary.files_partial} partial, {summary.files_failed} failed, "
                f"{summary.files_invalid} invalid, {summary.files_empty} empty "
                f"(see {os.path.join(output_dir, 'extraction_report.json')})"
            )
            for r in problems:
                logger.warning(
                    f"  {os.path.basename(r.file_path)}: {r.status.value}"
                    f"{' - ' + r.error_message if r.error_message else ''}"
                )
        success = n_ok > 0 and os.path.exists(product)
        msg = ""
        if not success:
            msg = f"no file produced a usable table (see {os.path.join(output_dir, 'extraction_report.json')})"
            logger.error(f"Extract: {msg}")
        elif not problems:
            logger.info("All processing complete!")
        return finish(product, rows, success, msg)
    except Exception as e:
        logger.exception(f"Extract failed: {e}")
        return finish(error=f"{type(e).__name__}: {e}")


# ===============================
# Main function
# ===============================
def main(argv=None):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        encoding="utf-8",
    )
    parser = argparse.ArgumentParser(
        description="GenBank Metadata Extractor (refactored): 3-strategy Assembly + record tracking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, nargs="+", help="Input GenBank files or directories")
    parser.add_argument("-o", "--output", default=".", help="Output directory")
    parser.add_argument("-m", "--metadata-output", default="metadata.csv", help="Metadata table (one row per record)")
    parser.add_argument("-a", "--assembly-output", default="assembly.csv", help="Assembly / sequencing data table")
    parser.add_argument("-f", "--final-output", default="final.csv", help="Metadata and assembly data merged")

    parser.add_argument(
        "--stream",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Stream processing (default; same results, less memory on large inputs)",
    )
    parser.add_argument("--batch", action="store_true", help="Batch mode (parallel)")
    parser.add_argument("--max_tasks", type=int, default=4, help="Parallel processes in --batch mode")
    parser.add_argument("--log_interval", type=int, default=5000, help="Log progress every N records (0 = never)")

    parser.add_argument("--metadata_delimiter", default=",", help="Column delimiter of the metadata table")
    parser.add_argument("--assembly_delimiter", default="\t", help="Column delimiter of the assembly table")
    parser.add_argument("--final_delimiter", default=",", help="Column delimiter of the merged table")

    exclusive_group = parser.add_mutually_exclusive_group()
    exclusive_group.add_argument(
        "--only_assembly_meta", "-oam", action="store_true", help="Extract only the assembly / sequencing data"
    )
    exclusive_group.add_argument(
        "--only_except_assembly_meta", "-oeam", action="store_true", help="Extract only the metadata (no assembly data)"
    )

    parser.add_argument(
        "--update_taxonomy", action="store_true", help="Update the local ete3 NCBI taxonomy database first"
    )
    parser.add_argument(
        "--no_report", action="store_true", help="Do not write extraction_report.json and the error table"
    )
    parser.add_argument("--no_tracking", action="store_true", help="Do not write the per-record tracking table")
    parser.add_argument("--no_recursive", action="store_true", help="Do not search subdirectories of the inputs")
    parser.add_argument("--if_use_chunksize", action="store_true", help="(Legacy) same as --stream")

    args = parser.parse_args(argv)
    result = extract(
        args.input,
        args.output,
        stream=args.stream or args.if_use_chunksize,
        batch=args.batch,
        max_tasks=args.max_tasks,
        final_name=args.final_output,
        metadata_name=args.metadata_output,
        assembly_name=args.assembly_output,
        metadata_delimiter=args.metadata_delimiter,
        assembly_delimiter=args.assembly_delimiter,
        final_delimiter=args.final_delimiter,
        only_assembly=args.only_assembly_meta,
        only_metadata=args.only_except_assembly_meta,
        recursive=not args.no_recursive,
        log_interval=args.log_interval,
        update_taxonomy=args.update_taxonomy,
        report=not args.no_report,
        details=True,
        tracking=not args.no_tracking,
    )
    if not result.success:
        sys.exit(1)


if __name__ == "__main__":
    main()
