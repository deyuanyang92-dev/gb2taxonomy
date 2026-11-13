"""Reading files record by record (stream or memory) and writing the tables."""

from __future__ import annotations

import csv
import logging
import os
import time
from typing import Any, Callable, Optional

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from g2t._extract.files import format_size, get_memory_mb, open_genbank
from g2t._extract.parse import (
    ASSEMBLY_CSV_COLUMNS,
    canonical_metadata_columns,
    extract_assembly_from_record,
    extract_metadata,
    prescan_fieldnames,
)
from g2t._extract.status import FileResult, FileStatus, ProcessingSummary, RecordError, RecordTrack
from g2t._extract.taxonomy import prefetch_taxonomy

logger = logging.getLogger("g2t.extract")  # one logger name for the whole step


# ===============================
# Unified record processing
# ===============================
def process_record_unified(
    record: SeqRecord,
    include_taxonomy: bool = True,
    do_metadata: bool = True,
    do_assembly: bool = True,
) -> tuple[dict[str, str] | None, dict[str, str] | None, str, bool]:
    metadata = None
    assembly_data = None
    assembly_source = ""
    has_asm = False

    if do_metadata:
        metadata = extract_metadata(record, include_taxonomy)

    if do_assembly:
        assembly_data, assembly_source = extract_assembly_from_record(record)
        has_asm = assembly_data is not None
        if assembly_data is not None and "ACCESSION" not in assembly_data:
            accs: Any = record.annotations.get("accessions") or [""]
            acc = str(accs[0] or record.id)
            assembly_data["ACCESSION"] = acc

    return metadata, assembly_data, assembly_source, has_asm


# ===============================
# Shared record processing loop
# ===============================
def count_record_terminators(file_path: str) -> int | None:
    """Number of lines that are exactly "//" (GenBank record terminators); None if unreadable."""
    count = 0
    try:
        with open_genbank(file_path, binary=True) as f:
            for line in f:
                if line.strip() == b"//":
                    count += 1
    except OSError:
        return None
    return count


RecordCallback = Callable[
    [Optional[dict[str, str]], Optional[dict[str, str]], str, bool, RecordTrack, "FileResult"], None
]


def _process_record_loop(
    file_path: str,
    callback: RecordCallback,
    include_taxonomy: bool = True,
    do_metadata: bool = True,
    do_assembly: bool = True,
    log_interval: int = 5000,
) -> FileResult:
    """Shared record processing loop for a single GB file.

    callback receives (metadata, assembly_data, assembly_source, has_asm, track, result) for each record.

    Optimized: Single-pass file reading (no pre-validation).
    """
    result = FileResult(file_path=file_path)
    start = time.monotonic()

    # Quick file existence check
    if not os.path.exists(file_path):
        result.status = FileStatus.INVALID
        result.error_message = "File does not exist"
        result.elapsed_time = time.monotonic() - start
        return result

    sz = os.path.getsize(file_path)
    if sz == 0:
        result.status = FileStatus.EMPTY
        result.error_message = "Empty file (0 bytes)"
        result.elapsed_time = time.monotonic() - start
        return result

    logger.info(f"  Size: {format_size(sz)}")

    record_idx = 0
    parse_error = ""
    try:
        with open_genbank(file_path) as handle:
            for record in SeqIO.parse(handle, "genbank"):
                record_idx += 1
                accession = ""
                track = RecordTrack(accession="", record_index=record_idx)

                try:
                    accession = record.annotations.get("accessions", [""])[0] or record.id
                    track.accession = accession

                    metadata, assembly_data, assembly_source, has_asm = process_record_unified(
                        record, include_taxonomy, do_metadata, do_assembly
                    )
                    if do_assembly:
                        # must be set before the callback so that a failing
                        # callback is counted as assembly_failed instead of vanishing.
                        track.has_assembly_data = has_asm

                    callback(metadata, assembly_data, assembly_source, has_asm, track, result)

                    if do_metadata and metadata is not None:
                        result.metadata_extracted += 1
                        track.metadata_ok = True
                    if do_assembly:
                        if assembly_data is not None:
                            result.assembly_extracted += 1
                            track.assembly_ok = True
                            track.assembly_source = assembly_source
                        elif not has_asm:
                            result.assembly_no_data += 1

                except Exception as e:
                    track.error_message = str(e)
                    if do_metadata and not track.metadata_ok:
                        result.metadata_failed += 1
                        result.add_record_error(
                            RecordError(
                                file_path=file_path,
                                record_index=record_idx,
                                accession=accession,
                                error_type=type(e).__name__,
                                error_message=str(e),
                                phase="metadata",
                            )
                        )
                    if do_assembly and track.has_assembly_data and not track.assembly_ok:
                        result.assembly_failed += 1

                result.record_tracks.append(track)

                # Progress logging
                done = result.metadata_extracted + result.metadata_failed
                if log_interval > 0 and done > 0 and done % log_interval == 0:
                    elapsed = time.monotonic() - start
                    speed = done / elapsed if elapsed > 0 else 0
                    mem = get_memory_mb()
                    mem_s = f", mem:{mem:.0f}MB" if mem else ""
                    logger.info(f"  Progress: {done} ({speed:.0f}/s{mem_s})")

    except Exception as e:
        parse_error = f"Biopython parse error: {e}"
        result.error_message = parse_error
        logger.error(f"  File-level error: {e}")
        result.add_record_error(
            RecordError(
                file_path=file_path,
                record_index=record_idx + 1,
                accession="",
                error_type=type(e).__name__,
                error_message=str(e),
                phase="parse",
            )
        )

    # Set total records after processing
    result.total_records = record_idx

    # reconcile the number of parsed records with the number of "//"
    # terminators, so records lost after a mid-file parse error (or a truncated
    # download) can never go unnoticed.
    file_problem = bool(parse_error)
    expected = count_record_terminators(file_path)
    result.records_expected = expected
    if expected is not None and expected != record_idx:
        file_problem = True
        name = os.path.basename(file_path)
        if record_idx > expected:
            detail = "the last record has no '//' terminator (file truncated?)"
        else:
            detail = f"{expected - record_idx} record(s) were not parsed"
        msg = (
            f"record count mismatch: file has {expected} '//' terminator(s) "
            f"but {record_idx} record(s) were parsed; {detail}"
        )
        logger.warning(f"  {name}: {msg}")
        result.error_message = f"{result.error_message}; {msg}" if result.error_message else msg
        result.add_record_error(
            RecordError(
                file_path=file_path,
                record_index=record_idx + 1,
                accession="",
                error_type="RecordCountMismatch",
                error_message=msg,
                phase="parse",
            )
        )

    # Status determination (a file with any parse problem is never SUCCESS)
    if do_metadata:
        if result.metadata_extracted == 0:
            result.status = FileStatus.FAILED
            if not result.error_message:
                result.error_message = "No metadata extracted"
        elif result.metadata_failed > 0 or file_problem:
            result.status = FileStatus.PARTIAL
        else:
            result.status = FileStatus.SUCCESS
    elif file_problem:
        done_ok = result.assembly_extracted + result.assembly_no_data
        result.status = FileStatus.PARTIAL if done_ok > 0 else FileStatus.FAILED
    else:
        result.status = FileStatus.SUCCESS

    result.elapsed_time = time.monotonic() - start
    return result


# ===============================
# Stream mode wrapper
# ===============================
def stream_process_single_file(
    file_path: str,
    metadata_csv_path: str | None,
    assembly_csv_path: str | None,
    include_taxonomy: bool = True,
    do_extract_metadata: bool = True,
    do_extract_assembly: bool = True,
    metadata_delimiter: str = ",",
    assembly_delimiter: str = "\t",
    log_interval: int = 5000,
) -> FileResult:
    metadata_fh = None
    metadata_writer = None
    metadata_header_written = False
    assembly_fh = None
    assembly_writer = None

    try:
        if do_extract_metadata and metadata_csv_path:
            metadata_fh = open(metadata_csv_path, "w", newline="", encoding="utf-8")  # noqa: SIM115
        if do_extract_assembly and assembly_csv_path:
            assembly_fh = open(assembly_csv_path, "w", newline="", encoding="utf-8")  # noqa: SIM115
            assembly_writer = csv.DictWriter(
                assembly_fh, fieldnames=ASSEMBLY_CSV_COLUMNS, delimiter=assembly_delimiter, extrasaction="ignore"
            )
            assembly_writer.writeheader()

        # Pre-scan to discover all source qualifier keys
        all_fieldnames = prescan_fieldnames([file_path]) if do_extract_metadata else None

        def stream_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
            nonlocal metadata_header_written, metadata_writer
            # Write metadata
            if do_extract_metadata and metadata is not None and metadata_fh:
                if not metadata_header_written:
                    metadata_writer = csv.DictWriter(
                        metadata_fh, fieldnames=all_fieldnames, delimiter=metadata_delimiter, extrasaction="ignore"
                    )
                    metadata_writer.writeheader()
                    metadata_header_written = True
                for fn in metadata_writer.fieldnames:
                    if fn not in metadata:
                        metadata[fn] = ""
                metadata_writer.writerow(metadata)

            # Write assembly
            if do_extract_assembly and assembly_data is not None and assembly_fh:
                for col in ASSEMBLY_CSV_COLUMNS:
                    if col not in assembly_data:
                        assembly_data[col] = ""
                assembly_writer.writerow(assembly_data)

        result = _process_record_loop(
            file_path, stream_callback, include_taxonomy, do_extract_metadata, do_extract_assembly, log_interval
        )
        return result

    finally:
        for fh in (metadata_fh, assembly_fh):
            if fh:
                try:
                    fh.flush()
                    fh.close()
                except Exception:
                    pass


# ===============================
# Memory mode wrapper
# ===============================
def memory_process_single_file(
    file_path: str,
    include_taxonomy: bool = True,
    do_extract_assembly: bool = True,
    do_extract_metadata: bool = True,
) -> tuple[FileResult, list[dict], list[dict]]:
    md_rows = []
    asm_rows = []

    def memory_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
        if do_extract_metadata and metadata:
            md_rows.append(metadata)
        if do_extract_assembly and assembly_data:
            asm_rows.append(assembly_data)

    result = _process_record_loop(
        file_path, memory_callback, include_taxonomy, do_extract_metadata, do_extract_assembly
    )
    return result, md_rows, asm_rows


# ===============================
# Merge metadata + assembly
# ===============================
def combine_and_save_final_csv(
    md_file: str, asm_file: str, out_file: str, delimiter: str = ",", md_delimiter: str = ",", asm_delimiter: str = "\t"
):
    start = time.monotonic()
    # read everything as text so leading zeros ("007"), exponent-like ids
    # ("1E5") and big integers are never rewritten; empty cells stay "" and
    # strings such as "NA" are not turned into NaN.
    try:
        md_df = pd.read_csv(md_file, sep=md_delimiter, dtype=str, keep_default_na=False, low_memory=False)
    except Exception as e:
        logger.error(f"Failed to read metadata: {e}")
        return

    asm_df = None
    try:
        asm_df = pd.read_csv(asm_file, sep=asm_delimiter, dtype=str, keep_default_na=False)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        pass
    except Exception as e:
        logger.warning(f"Failed to read assembly: {e}")

    if asm_df is not None and "ACCESSION" in md_df.columns and "ACCESSION" in asm_df.columns:
        md_df["ACCESSION"] = md_df["ACCESSION"].astype(str).str.strip()
        asm_df["ACCESSION"] = asm_df["ACCESSION"].astype(str).str.strip()
        asm_df = asm_df.drop_duplicates(subset="ACCESSION", keep="first")
        final = pd.merge(md_df, asm_df, how="left", on="ACCESSION").fillna("")
    else:
        final = md_df

    final.to_csv(out_file, index=False, sep=delimiter, encoding="utf-8")
    logger.info(f"Merged: {out_file} ({len(final)} records, {time.monotonic() - start:.1f}s)")


def _write_metadata_frame(rows: list[dict], path: str, delimiter: str) -> pd.DataFrame:
    """Write metadata rows with the same column set/order as the stream writer."""
    df = pd.DataFrame(rows)
    df = df.reindex(columns=canonical_metadata_columns(df.columns)).fillna("")
    df.to_csv(path, index=False, sep=delimiter, encoding="utf-8")
    return df


def _write_assembly_frame(rows: list[dict], path: str, delimiter: str) -> pd.DataFrame:
    """Write assembly rows (header-only when empty, like the stream writer)."""
    df = pd.DataFrame(rows)
    for c in ASSEMBLY_CSV_COLUMNS:
        if c not in df.columns:
            df[c] = ""
    df = df[ASSEMBLY_CSV_COLUMNS].fillna("")
    df.to_csv(path, index=False, sep=delimiter, encoding="utf-8")
    return df


# ===============================
# Main flow: non-batch
# ===============================
def process_all_files(
    gb_files: list[str],
    output_dir: str,
    metadata_output: str = "metadata.csv",
    assembly_output: str = "assembly.csv",
    final_output: str = "final.csv",
    include_taxonomy: bool = True,
    only_assembly: bool = False,
    only_except_assembly: bool = False,
    metadata_delimiter: str = ",",
    assembly_delimiter: str = "\t",
    final_delimiter: str = ",",
    stream_mode: bool = False,
    log_interval: int = 5000,
) -> ProcessingSummary:
    summary = ProcessingSummary()
    summary.gb_files_found = len(gb_files)
    t0 = time.monotonic()

    do_metadata = not only_assembly
    do_assembly = not only_except_assembly
    tax = include_taxonomy and not only_except_assembly

    metadata_path = os.path.join(output_dir, metadata_output)
    assembly_path = os.path.join(output_dir, assembly_output)
    final_path = os.path.join(output_dir, final_output)

    prefetch_taxonomy(gb_files, output_dir, tax)

    if stream_mode:
        metadata_fh = None
        metadata_writer = None
        metadata_header_written = False
        assembly_fh = None
        assembly_writer = None
        try:
            if do_metadata:
                metadata_fh = open(metadata_path, "w", newline="", encoding="utf-8")  # noqa: SIM115
            if do_assembly:
                assembly_fh = open(assembly_path, "w", newline="", encoding="utf-8")  # noqa: SIM115
                assembly_writer = csv.DictWriter(
                    assembly_fh, fieldnames=ASSEMBLY_CSV_COLUMNS, delimiter=assembly_delimiter, extrasaction="ignore"
                )
                assembly_writer.writeheader()

            # Pre-scan all files to discover every source qualifier key
            all_fieldnames = prescan_fieldnames(gb_files) if do_metadata else None

            for i, fp in enumerate(gb_files, 1):
                logger.info(f"[{i}/{len(gb_files)}] {os.path.basename(fp)} ({format_size(os.path.getsize(fp))})")

                def stream_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
                    nonlocal metadata_header_written, metadata_writer
                    if do_metadata and metadata and metadata_fh:
                        if not metadata_header_written:
                            metadata_writer = csv.DictWriter(
                                metadata_fh,
                                fieldnames=all_fieldnames,
                                delimiter=metadata_delimiter,
                                extrasaction="ignore",
                            )
                            metadata_writer.writeheader()
                            metadata_header_written = True
                        for fn in metadata_writer.fieldnames:
                            if fn not in metadata:
                                metadata[fn] = ""
                        metadata_writer.writerow(metadata)
                    if do_assembly and assembly_data and assembly_fh:
                        for c in ASSEMBLY_CSV_COLUMNS:
                            if c not in assembly_data:
                                assembly_data[c] = ""
                        assembly_writer.writerow(assembly_data)

                fr = _process_record_loop(fp, stream_callback, tax, do_metadata, do_assembly, log_interval)
                summary.add_result(fr)

                icon = {"success": "OK", "partial": "~", "failed": "X"}.get(fr.status.value, "?")
                logger.info(
                    f"  {icon} MD:{fr.metadata_extracted}/{fr.total_records} "
                    f"ASM:{fr.assembly_extracted} ({fr.elapsed_time:.1f}s)"
                )

        finally:
            for fh in (metadata_fh, assembly_fh):
                if fh:
                    try:
                        fh.flush()
                        fh.close()
                    except Exception:
                        pass

        if do_metadata and do_assembly and summary.total_metadata_extracted > 0:
            try:
                combine_and_save_final_csv(
                    metadata_path, assembly_path, final_path, final_delimiter, metadata_delimiter, assembly_delimiter
                )
            except Exception as e:
                logger.error(f"Merge error: {e}")

    else:
        all_metadata = []
        all_assembly = []
        for i, fp in enumerate(gb_files, 1):
            logger.info(f"[{i}/{len(gb_files)}] {os.path.basename(fp)}")

            def memory_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
                if do_metadata and metadata:
                    all_metadata.append(metadata)
                if do_assembly and assembly_data:
                    all_assembly.append(assembly_data)

            fr = _process_record_loop(fp, memory_callback, tax, do_metadata, do_assembly, log_interval)
            summary.add_result(fr)

            icon = {"success": "OK", "partial": "~", "failed": "X"}.get(fr.status.value, "?")
            logger.info(
                f"  {icon} MD:{fr.metadata_extracted}/{fr.total_records} "
                f"ASM:{fr.assembly_extracted} ({fr.elapsed_time:.1f}s)"
            )

        if do_metadata and all_metadata:
            try:
                _write_metadata_frame(all_metadata, metadata_path, metadata_delimiter)
                logger.info(f"Metadata: {metadata_path} ({len(all_metadata)} records)")
            except Exception as e:
                logger.error(f"Failed to save metadata: {e}")

        if do_assembly and (all_assembly or all_metadata):
            try:
                _write_assembly_frame(all_assembly, assembly_path, assembly_delimiter)
                logger.info(f"Assembly: {assembly_path} ({len(all_assembly)} records)")
            except Exception as e:
                logger.error(f"Failed to save assembly: {e}")

        if do_metadata and do_assembly and all_metadata:
            try:
                combine_and_save_final_csv(
                    metadata_path, assembly_path, final_path, final_delimiter, metadata_delimiter, assembly_delimiter
                )
            except Exception as e:
                logger.error(f"Merge error: {e}")

    summary.total_elapsed_time = time.monotonic() - t0
    return summary
