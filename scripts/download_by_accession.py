#!/usr/bin/env python3
"""
Download GenBank records by accession number list.

Usage:
    python scripts/download_by_accession.py -i accession_list.txt -e you@example.com -k API_KEY
    python scripts/download_by_accession.py -i accession_list.txt -e you@example.com -k API_KEY --resume
"""

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

from Bio import Entrez


def parse_args():
    p = argparse.ArgumentParser(description="Download GenBank records by accession list")
    p.add_argument("-i", "--input", required=True, help="File with one accession per line")
    p.add_argument("-o", "--output", default=None, help="Output directory (default: <input>_gb/)")
    p.add_argument("-b", "--batch-size", type=int, default=200, help="Accessions per batch (default: 200)")
    p.add_argument("-e", "--email", required=True, help="Email for NCBI")
    p.add_argument("-k", "--api-key", default=None, help="NCBI API key")
    p.add_argument("--resume", action="store_true", help="Resume from previous download")
    return p.parse_args()


def load_accessions(path: Path) -> list[str]:
    accs = []
    with open(path) as f:
        for line in f:
            acc = line.strip()
            if acc:
                accs.append(acc)
    return accs


def load_progress(output_dir: Path) -> dict:
    pfile = output_dir / "download_progress.json"
    if pfile.exists():
        with open(pfile) as f:
            return json.load(f)
    return {}


def save_progress(output_dir: Path, progress: dict):
    pfile = output_dir / "download_progress.json"
    progress["timestamp"] = datetime.now().isoformat()
    with open(pfile, "w") as f:
        json.dump(progress, f, indent=2)


def download_batch(accessions: list[str], output_path: Path) -> bool:
    try:
        handle = Entrez.efetch(db="nuccore", id=accessions, rettype="gb", retmode="text")
        text = handle.read()
        handle.close()
        if not text.strip():
            return False
        with open(output_path, "w") as f:
            f.write(text)
        return True
    except Exception as e:
        print(f"  Error: {e}", file=sys.stderr)
        return False


def main():
    args = parse_args()
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    input_path = Path(args.input)
    accessions = load_accessions(input_path)
    output_dir = Path(args.output) if args.output else input_path.parent / f"{input_path.stem}_gb"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Accessions: {len(accessions)}")
    print(f"Output: {output_dir.resolve()}")

    progress = load_progress(output_dir)
    if args.resume and progress:
        start_batch = progress.get("last_batch", -1) + 1
        print(f"Resuming from batch {start_batch + 1}")
    else:
        start_batch = 0

    num_batches = (len(accessions) + args.batch_size - 1) // args.batch_size
    delay = 0.1 if args.api_key else 0.34
    failed = 0

    for i in range(start_batch, num_batches):
        start = i * args.batch_size
        end = min(start + args.batch_size, len(accessions))
        batch = accessions[start:end]
        batch_file = output_dir / f"batch_{i + 1:04d}.gb"

        print(f"  Batch {i + 1}/{num_batches} ({start + 1}-{end}/{len(accessions)})...", end=" ", flush=True)

        if download_batch(batch, batch_file):
            size_kb = batch_file.stat().st_size / 1024
            print(f"OK ({size_kb:.0f} KB)")
        else:
            print("FAILED")
            failed += 1
            if failed >= 3:
                print("Too many failures. Use --resume to continue.")
                progress.update({"last_batch": i - 1, "completed": False})
                save_progress(output_dir, progress)
                sys.exit(1)
            continue

        progress.update({"last_batch": i, "completed": i == num_batches - 1})
        save_progress(output_dir, progress)
        time.sleep(delay)

    # Count LOCUS lines to verify
    total_locus = 0
    for f in sorted(output_dir.glob("batch_*.gb")):
        with open(f) as fh:
            total_locus += sum(1 for line in fh if line.startswith("LOCUS"))

    gb_files = list(output_dir.glob("batch_*.gb"))
    total_size = sum(f.stat().st_size for f in gb_files) / (1024 * 1024)
    print(f"\nDone! {len(gb_files)} files, {total_locus} records, {total_size:.1f} MB")


if __name__ == "__main__":
    main()
