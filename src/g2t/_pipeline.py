"""Pipeline orchestrator — calls step functions directly (no subprocess)."""

from __future__ import annotations

import importlib
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from g2t.utils import is_valid_output


@dataclass
class PipelineResult:
    success: bool
    output_dir: str
    steps_completed: int = 0
    final_output: str = ""
    elapsed: float = 0.0
    log: List[str] = field(default_factory=list)


def _run_step(
    *,
    step_num: int,
    skip: bool,
    resume: bool,
    output_path: Path,
    module_path: str,
    func_name: str,
    kwargs: dict,
    quiet: bool,
    log_lines: List[str],
    fallback_dir: Optional[Path] = None,
    fallback_glob: Optional[str] = None,
) -> Optional[str]:
    """Run a single pipeline step.

    Returns an error message string on failure, None on success.
    Raises FileNotFoundError if fallback check fails.
    """
    if skip:
        if fallback_dir and fallback_glob and not is_valid_output(output_path):
            _find_fallback(fallback_dir, fallback_glob, output_path, f"Step {step_num}")
        return None

    if resume and is_valid_output(output_path):
        msg = f"Step {step_num} skipped (resume): {output_path}"
        log_lines.append(msg)
        if not quiet:
            print(msg)
    else:
        mod = importlib.import_module(module_path)
        step_func = getattr(mod, func_name)
        result = step_func(**kwargs)
        if not result.success:
            detail = getattr(result, "output_file", "")
            return f"Step {step_num} failed" + (f": {detail}" if detail else "")
        msg = f"Step {step_num} done: {output_path} ({result.rows} rows, {result.elapsed:.1f}s)"
        log_lines.append(msg)
        if not quiet:
            print(msg)

    if fallback_dir and fallback_glob and not is_valid_output(output_path):
        _find_fallback(fallback_dir, fallback_glob, output_path, f"Step {step_num}")
    return None


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

    def _fail(msg: str) -> PipelineResult:
        return PipelineResult(
            success=False, output_dir=str(out_root),
            steps_completed=steps_completed, elapsed=time.time() - start_time,
            log=log_lines + [msg],
        )

    # --- Step 1: Extract ---
    step1_kwargs = dict(
        input_files=input_files,
        output_dir=str(extract_dir),
        stream=stream,
        batch=batch,
        max_tasks=max_tasks,
        final_name="final.csv",
    )
    if extract_extra:
        step1_kwargs.update(extract_extra)
    err = _run_step(
        step_num=1, skip=skip_extract, resume=resume,
        output_path=out_final_csv, module_path="g2t.extract", func_name="extract",
        kwargs=step1_kwargs, quiet=quiet, log_lines=log_lines,
        fallback_dir=extract_dir, fallback_glob="*final*.csv",
    )
    if err:
        return _fail(err)
    if not skip_extract:
        steps_completed += 1

    # --- Step 2: Classify ---
    step2_kwargs = dict(
        input_file=str(out_final_csv),
        output_dir=str(classify_dir),
    )
    if classify_extra:
        step2_kwargs.update(classify_extra)
    err = _run_step(
        step_num=2, skip=skip_classify, resume=resume,
        output_path=out_assigned_all, module_path="g2t.classify", func_name="classify",
        kwargs=step2_kwargs, quiet=quiet, log_lines=log_lines,
        fallback_dir=classify_dir, fallback_glob="*assigned*all*.csv",
    )
    if err:
        return _fail(err)
    if not skip_classify:
        steps_completed += 1

    # --- Step 3: Voucher ---
    step3_kwargs = dict(
        file_path=str(out_assigned_all),
        output_dir=str(voucher_dir),
    )
    if normalize_columns:
        step3_kwargs["normalize_column_names"] = True
    if voucher_extra:
        step3_kwargs.update(voucher_extra)
    err = _run_step(
        step_num=3, skip=skip_voucher, resume=resume,
        output_path=out_updated_voucher, module_path="g2t.voucher",
        func_name="build_species_vouchers",
        kwargs=step3_kwargs, quiet=quiet, log_lines=log_lines,
        fallback_dir=voucher_dir, fallback_glob="updated_species_voucher*.csv",
    )
    if err:
        return _fail(err)
    if not skip_voucher:
        steps_completed += 1

    # --- Step 4: Organize ---
    step4_kwargs = dict(
        input_file=str(out_updated_voucher),
        output_file=str(out_organized),
    )
    if organize_extra:
        step4_kwargs.update(organize_extra)
    err = _run_step(
        step_num=4, skip=skip_organize, resume=resume,
        output_path=out_organized, module_path="g2t.organize", func_name="organize",
        kwargs=step4_kwargs, quiet=quiet, log_lines=log_lines,
    )
    if err:
        return _fail(err)
    if not skip_organize:
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
