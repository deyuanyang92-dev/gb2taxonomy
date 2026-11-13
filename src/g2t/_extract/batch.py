"""Batch mode: one process per file, then merge the per-file tables."""

from __future__ import annotations

import concurrent.futures
import functools
import logging
import multiprocessing
import os
import time

import pandas as pd

from g2t._extract import taxonomy
from g2t._extract.parse import canonical_metadata_columns
from g2t._extract.process import (
    _write_assembly_frame,
    _write_metadata_frame,
    combine_and_save_final_csv,
    memory_process_single_file,
    stream_process_single_file,
)
from g2t._extract.status import FileResult, FileStatus, ProcessingSummary
from g2t._extract.taxonomy import _init_worker, prefetch_taxonomy

logger = logging.getLogger("g2t.extract")  # one logger name for the whole step


# ===============================
# Batch mode
# ===============================
def _unique_batch_names(gb_files: list[str]) -> list[str]:
    """One output-directory stem per input file; same-named files get a numeric suffix."""
    used = set()
    names = []
    for fp in gb_files:
        stem = os.path.splitext(os.path.basename(fp))[0]
        name, k = stem, 1
        while name.lower() in used:
            k += 1
            name = f"{stem}_{k}"
        used.add(name.lower())
        names.append(name)
    return names


def _batch_paths(output_dir: str, base: str) -> tuple[str, str, str, str]:
    sub = os.path.join(output_dir, base + "_metadata")
    return (
        sub,
        os.path.join(sub, base + "_metadata.csv"),
        os.path.join(sub, base + "_assembly.csv"),
        os.path.join(sub, base + "_metadata_final.csv"),
    )


def _process_single_file_batch(
    file_path: str,
    base: str,
    *,
    output_dir: str,
    do_metadata: bool,
    do_assembly: bool,
    tax: bool,
    stream_mode: bool,
    metadata_delimiter: str,
    assembly_delimiter: str,
    final_delimiter: str,
    log_interval: int,
) -> FileResult:
    """Process one GenBank file into ``<output_dir>/<base>_metadata/``.

    Module-level on purpose: ProcessPoolExecutor must be able to pickle it;
    options are bound with ``functools.partial``.
    """
    sub, md_name, asm_name, final_path = _batch_paths(output_dir, base)
    os.makedirs(sub, exist_ok=True)
    metadata_path = md_name if do_metadata else None
    assembly_path = asm_name if do_assembly else None

    if stream_mode:
        fr = stream_process_single_file(
            file_path,
            metadata_path,
            assembly_path,
            tax,
            do_metadata,
            do_assembly,
            metadata_delimiter,
            assembly_delimiter,
            log_interval,
        )
    else:
        fr, mdr, asmr = memory_process_single_file(file_path, tax, do_assembly, do_metadata)
        if do_metadata and mdr and metadata_path:
            try:
                _write_metadata_frame(mdr, metadata_path, metadata_delimiter)
            except Exception as e:
                logger.error(f"Save error: {e}")
        if do_assembly and asmr and assembly_path:
            try:
                _write_assembly_frame(asmr, assembly_path, assembly_delimiter)
            except Exception as e:
                logger.error(f"Save error: {e}")

    if (
        do_metadata
        and do_assembly
        and metadata_path
        and assembly_path
        and os.path.exists(metadata_path)
        and os.path.getsize(metadata_path) > 0
    ):
        try:
            combine_and_save_final_csv(
                metadata_path, assembly_path, final_path, final_delimiter, metadata_delimiter, assembly_delimiter
            )
        except Exception as e:
            logger.error(f"Merge error: {e}")

    return fr


def _read_text_csv(path: str, sep: str) -> pd.DataFrame | None:
    """Read a CSV as text; None when missing/empty."""
    try:
        return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, low_memory=False)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None
    except Exception as e:
        logger.warning(f"Failed to read {path}: {e}")
        return None


def _merge_batch_outputs(
    names: list[str],
    output_dir: str,
    metadata_output: str,
    assembly_output: str,
    final_output: str,
    do_metadata: bool,
    do_assembly: bool,
    metadata_delimiter: str,
    assembly_delimiter: str,
    final_delimiter: str,
) -> int:
    """Concatenate the per-file CSVs (in input order) into metadata/assembly/final.csv.

    The result has the same columns and order as the non-batch output.
    Returns the number of merged metadata rows.
    """
    metadata_path = os.path.join(output_dir, metadata_output)
    assembly_path = os.path.join(output_dir, assembly_output)
    final_path = os.path.join(output_dir, final_output)

    md_frames, asm_frames = [], []
    for base in names:
        _, md_file, asm_file, _ = _batch_paths(output_dir, base)
        if do_metadata:
            df = _read_text_csv(md_file, metadata_delimiter)
            if df is not None and len(df):
                md_frames.append(df)
        if do_assembly:
            df = _read_text_csv(asm_file, assembly_delimiter)
            if df is not None and len(df):
                asm_frames.append(df)

    n_rows = 0
    if do_metadata and md_frames:
        md = pd.concat(md_frames, ignore_index=True, sort=False)
        md = md.reindex(columns=canonical_metadata_columns(md.columns)).fillna("")
        md.to_csv(metadata_path, index=False, sep=metadata_delimiter, encoding="utf-8")
        n_rows = len(md)
        logger.info(f"Metadata: {metadata_path} ({n_rows} records)")
    if do_assembly and (asm_frames or (md_frames and do_metadata)):
        rows = pd.concat(asm_frames, ignore_index=True, sort=False).to_dict("records") if asm_frames else []
        _write_assembly_frame(rows, assembly_path, assembly_delimiter)
        logger.info(f"Assembly: {assembly_path} ({len(rows)} records)")
    if do_metadata and do_assembly and n_rows > 0:
        combine_and_save_final_csv(
            metadata_path, assembly_path, final_path, final_delimiter, metadata_delimiter, assembly_delimiter
        )
    return n_rows


def process_batch_mode(
    gb_files: list[str],
    output_dir: str,
    max_tasks: int = 4,
    stream_mode: bool = False,
    only_assembly: bool = False,
    only_except_assembly: bool = False,
    include_taxonomy: bool = True,
    metadata_delimiter: str = ",",
    assembly_delimiter: str = "\t",
    final_delimiter: str = ",",
    log_interval: int = 5000,
    metadata_output: str = "metadata.csv",
    assembly_output: str = "assembly.csv",
    final_output: str = "final.csv",
) -> ProcessingSummary:
    summary = ProcessingSummary()
    summary.gb_files_found = len(gb_files)
    t0 = time.monotonic()
    os.makedirs(output_dir, exist_ok=True)

    do_metadata = not only_assembly
    do_assembly = not only_except_assembly
    tax = include_taxonomy and not only_except_assembly

    names = _unique_batch_names(gb_files)
    worker = functools.partial(
        _process_single_file_batch,
        output_dir=output_dir,
        do_metadata=do_metadata,
        do_assembly=do_assembly,
        tax=tax,
        stream_mode=stream_mode,
        metadata_delimiter=metadata_delimiter,
        assembly_delimiter=assembly_delimiter,
        final_delimiter=final_delimiter,
        log_interval=log_interval,
    )

    # bulk-load lineages once in the parent; workers only read the cache.
    prefetch_taxonomy(gb_files, output_dir, tax)

    # Shared taxonomy cache via multiprocessing.Manager
    manager = multiprocessing.Manager()
    shared_cache = manager.dict()
    # Pre-populate with current cache
    shared_cache.update(taxonomy._taxonomy.cache)

    logger.info(f"Batch mode: {len(gb_files)} files, concurrency {max_tasks}, stream={'yes' if stream_mode else 'no'}")

    results: dict[int, FileResult] = {}
    try:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max(1, max_tasks),
            initializer=_init_worker,
            initargs=(shared_cache, taxonomy._taxonomy.ete3_available, taxonomy._taxonomy.entrez_fallback),
        ) as ex:
            future_to_idx = {ex.submit(worker, f, names[i]): i for i, f in enumerate(gb_files)}
            for done, fut in enumerate(concurrent.futures.as_completed(future_to_idx), 1):
                idx = future_to_idx[fut]
                f = gb_files[idx]
                try:
                    fr = fut.result()
                    icon = {"success": "OK", "partial": "~", "failed": "X"}.get(fr.status.value, "?")
                    logger.info(
                        f"[{done}/{len(gb_files)}] {icon} {os.path.basename(f)} "
                        f"MD:{fr.metadata_extracted}/{fr.total_records} "
                        f"ASM:{fr.assembly_extracted} ({fr.elapsed_time:.1f}s)"
                    )
                except Exception as e:
                    fr = FileResult(file_path=f, status=FileStatus.FAILED, error_message=str(e))
                    logger.error(f"[{done}/{len(gb_files)}] X {os.path.basename(f)}: {e}")
                results[idx] = fr
        # Sync cache back
        taxonomy._taxonomy.cache.update(dict(shared_cache))
    finally:
        manager.shutdown()

    # Report in input order (deterministic), not completion order.
    for i in range(len(gb_files)):
        if i in results:
            summary.add_result(results[i])

    # merge per-file outputs into <output_dir>/final.csv (+ metadata.csv, assembly.csv)
    merged_names = [
        names[i]
        for i in range(len(gb_files))
        if i in results and results[i].status in (FileStatus.SUCCESS, FileStatus.PARTIAL)
    ]
    try:
        _merge_batch_outputs(
            merged_names,
            output_dir,
            metadata_output,
            assembly_output,
            final_output,
            do_metadata,
            do_assembly,
            metadata_delimiter,
            assembly_delimiter,
            final_delimiter,
        )
    except Exception as e:
        logger.error(f"Batch merge error: {e}")

    summary.total_elapsed_time = time.monotonic() - t0
    return summary
