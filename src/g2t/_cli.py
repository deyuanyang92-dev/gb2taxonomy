"""
CLI entry points for g2t.

Installed as console scripts:
  g2t             → main_pipeline
  g2t-extract     → main_extract
  g2t-classify    → main_classify
  g2t-voucher     → main_voucher
  g2t-organize    → main_organize
"""

import argparse
import sys


def main_pipeline(argv=None):
    """g2t — Run the full GenBank-to-Taxonomy pipeline."""
    ap = argparse.ArgumentParser(
        prog="g2t",
        description="GenBank to Taxonomy pipeline: extract -> classify -> voucher -> organize",
    )
    ap.add_argument("-i", "--input", nargs="+", required=True,
                    help="GenBank files or directories")
    ap.add_argument("-o", "--output", required=True, help="Root output directory")
    ap.add_argument("--stream", action="store_true",
                    help="Use streaming mode for extraction (recommended for large files)")
    ap.add_argument("--batch", action="store_true", help="Use batch parallel processing")
    ap.add_argument("--max_tasks", type=int, default=1, help="Max parallel tasks for batch mode")
    ap.add_argument("--resume", action="store_true",
                    help="Skip steps whose output already exists")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--skip_extract", action="store_true")
    ap.add_argument("--skip_classify", action="store_true")
    ap.add_argument("--skip_voucher", action="store_true")
    ap.add_argument("--skip_organize", action="store_true")
    ap.add_argument("--normalize_columns", action="store_true",
                    help="Normalize column names (spaces to underscores)")
    ap.add_argument("--extract_extra", default="", help="Extra args for extract step")
    ap.add_argument("--classify_extra", default="", help="Extra args for classify step")
    ap.add_argument("--voucher_extra", default="", help="Extra args for voucher step")
    ap.add_argument("--organize_extra", default="", help="Extra args for organize step")

    args = ap.parse_args(argv)

    from g2t._pipeline import run
    result = run(
        input_files=args.input,
        output_dir=args.output,
        stream=args.stream,
        batch=args.batch,
        max_tasks=args.max_tasks,
        resume=args.resume,
        quiet=args.quiet,
        skip_extract=args.skip_extract,
        skip_classify=args.skip_classify,
        skip_voucher=args.skip_voucher,
        skip_organize=args.skip_organize,
        normalize_columns=args.normalize_columns,
    )

    if not result.success:
        for line in result.log:
            print(line, file=sys.stderr)
        raise SystemExit(1)

    if not args.quiet:
        print(f"\nAll steps complete!")
        print(f"  Output: {result.output_dir}")
        print(f"  Steps:  {result.steps_completed}")
        print(f"  Time:   {result.elapsed:.1f}s")

    return result


def main_extract(argv=None):
    """g2t-extract — Extract metadata from GenBank files."""
    from g2t.extract import main
    return main(argv)


def main_classify(argv=None):
    """g2t-classify — Classify gene types from metadata CSV."""
    from g2t.classify import main
    return main(argv)


def main_voucher(argv=None):
    """g2t-voucher — Build species voucher identifiers."""
    from g2t.voucher import main
    return main(argv)


def main_organize(argv=None):
    """g2t-organize — Organize genes by species."""
    from g2t.organize import main
    return main(argv)
