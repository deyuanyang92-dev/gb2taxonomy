"""
Pipeline orchestrator — calls step functions directly (no subprocess).
"""

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from g2t.utils import is_valid_output


@dataclass
class PipelineResult:
    success: bool
    output_dir: str
    steps_completed: int = 0
    final_output: str = ""
    elapsed: float = 0.0
    log: List[str] = field(default_factory=list)


def run(
    input_files: List[str],
    output_dir: str,
    stream: bool = True,
    batch: bool = False,
    max_tasks: int = 1,
    resume: bool = False,
    quiet: bool = False,
    skip_extract: bool = False,
    skip_classify: bool = False,
    skip_voucher: bool = False,
    skip_organize: bool = False,
    normalize_columns: bool = False,
    extract_extra: dict = None,
    classify_extra: dict = None,
    voucher_extra: dict = None,
    organize_extra: dict = None,
) -> PipelineResult:
    """Run the full g2t pipeline: extract -> classify -> voucher -> organize.

    Returns a PipelineResult with status and output paths.
    """
    start_time = time.time()
    log_lines: List[str] = []
    steps_completed = 0

    out_root = Path(output_dir).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    extract_dir = out_root / "gb_metadata"
    classify_dir = out_root / "labeled_genes"
    voucher_dir = out_root / "updated_species_vouchers"
    organize_dir = out_root / "organized_genes"

    for d in (extract_dir, classify_dir, voucher_dir, organize_dir):
        d.mkdir(parents=True, exist_ok=True)

    out_final_csv = extract_dir / "final.csv"
    out_assigned_all = classify_dir / "assigned_genes_types_all.csv"
    out_updated_voucher = voucher_dir / "updated_species_voucher.csv"
    out_organized = organize_dir / "organized_species_voucher.csv"

    # --- Step 1: Extract ---
    if not skip_extract:
        if resume and is_valid_output(out_final_csv):
            msg = f"Step 1 skipped (resume): {out_final_csv}"
            log_lines.append(msg)
            if not quiet:
                print(msg)
        else:
            from g2t.extract import extract as step_extract

            step_kwargs = dict(
                input_files=input_files,
                output_dir=str(extract_dir),
                stream=stream,
                batch=batch,
                max_tasks=max_tasks,
                final_name="final.csv",
            )
            if extract_extra:
                step_kwargs.update(extract_extra)

            result = step_extract(**step_kwargs)
            if not result.success:
                return PipelineResult(
                    success=False, output_dir=str(out_root),
                    steps_completed=steps_completed, elapsed=time.time() - start_time,
                    log=log_lines + [f"Step 1 failed: {result.output_file}"],
                )
            msg = f"Step 1 done: {out_final_csv} ({result.rows} rows, {result.elapsed:.1f}s)"
            log_lines.append(msg)
            if not quiet:
                print(msg)

        if not is_valid_output(out_final_csv):
            _find_fallback(extract_dir, "*final*.csv", out_final_csv, "Step 1")
        steps_completed += 1
    else:
        if not is_valid_output(out_final_csv):
            _find_fallback(extract_dir, "*final*.csv", out_final_csv, "Step 1")

    # --- Step 2: Classify ---
    if not skip_classify:
        if resume and is_valid_output(out_assigned_all):
            msg = f"Step 2 skipped (resume): {out_assigned_all}"
            log_lines.append(msg)
            if not quiet:
                print(msg)
        else:
            from g2t.classify import classify as step_classify

            step_kwargs = dict(
                input_file=str(out_final_csv),
                output_dir=str(classify_dir),
            )
            if classify_extra:
                step_kwargs.update(classify_extra)

            result = step_classify(**step_kwargs)
            if not result.success:
                return PipelineResult(
                    success=False, output_dir=str(out_root),
                    steps_completed=steps_completed, elapsed=time.time() - start_time,
                    log=log_lines + [f"Step 2 failed"],
                )
            msg = f"Step 2 done: {out_assigned_all} ({result.rows} rows, {result.elapsed:.1f}s)"
            log_lines.append(msg)
            if not quiet:
                print(msg)

        if not is_valid_output(out_assigned_all):
            _find_fallback(classify_dir, "*assigned*all*.csv", out_assigned_all, "Step 2")
        steps_completed += 1
    else:
        if not is_valid_output(out_assigned_all):
            _find_fallback(classify_dir, "*assigned*all*.csv", out_assigned_all, "Step 2")

    # --- Step 3: Voucher ---
    if not skip_voucher:
        if resume and is_valid_output(out_updated_voucher):
            msg = f"Step 3 skipped (resume): {out_updated_voucher}"
            log_lines.append(msg)
            if not quiet:
                print(msg)
        else:
            from g2t.voucher import build_species_vouchers as step_voucher

            step_kwargs = dict(
                file_path=str(out_assigned_all),
                output_dir=str(voucher_dir),
            )
            if normalize_columns:
                step_kwargs["normalize_column_names"] = True
            if voucher_extra:
                step_kwargs.update(voucher_extra)

            result = step_voucher(**step_kwargs)
            if not result.success:
                return PipelineResult(
                    success=False, output_dir=str(out_root),
                    steps_completed=steps_completed, elapsed=time.time() - start_time,
                    log=log_lines + [f"Step 3 failed"],
                )
            msg = f"Step 3 done: {out_updated_voucher} ({result.rows} rows, {result.elapsed:.1f}s)"
            log_lines.append(msg)
            if not quiet:
                print(msg)

        if not is_valid_output(out_updated_voucher):
            _find_fallback(voucher_dir, "updated_species_voucher*.csv", out_updated_voucher, "Step 3")
        steps_completed += 1
    else:
        if not is_valid_output(out_updated_voucher):
            _find_fallback(voucher_dir, "updated_species_voucher*.csv", out_updated_voucher, "Step 3")

    # --- Step 4: Organize ---
    if not skip_organize:
        if resume and is_valid_output(out_organized):
            msg = f"Step 4 skipped (resume): {out_organized}"
            log_lines.append(msg)
            if not quiet:
                print(msg)
        else:
            from g2t.organize import organize as step_organize

            step_kwargs = dict(
                input_file=str(out_updated_voucher),
                output_file=str(out_organized),
            )
            if organize_extra:
                step_kwargs.update(organize_extra)

            result = step_organize(**step_kwargs)
            if not result.success:
                return PipelineResult(
                    success=False, output_dir=str(out_root),
                    steps_completed=steps_completed, elapsed=time.time() - start_time,
                    log=log_lines + [f"Step 4 failed"],
                )
            msg = f"Step 4 done: {out_organized} ({result.rows} rows, {result.elapsed:.1f}s)"
            log_lines.append(msg)
            if not quiet:
                print(msg)

        steps_completed += 1

    elapsed = time.time() - start_time
    return PipelineResult(
        success=True,
        output_dir=str(out_root),
        steps_completed=steps_completed,
        final_output=str(out_organized),
        elapsed=elapsed,
        log=log_lines,
    )


def _find_fallback(directory: Path, pattern: str, target: Path, step_name: str):
    """If the expected output doesn't exist, search for a fallback."""
    if is_valid_output(target):
        return
    candidates = sorted(directory.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
    if candidates and candidates[0].stat().st_size > 0:
        return
    raise FileNotFoundError(
        f"{step_name}: no valid output found matching {pattern} in {directory}"
    )
