#!/usr/bin/env python3
"""
Download GenBank records by taxon using Biopython Entrez.

Usage:
    # By taxon name (any rank: order, family, genus, species...)
    python scripts/download_entrez.py --taxon Phyllodocida --email you@example.com --api-key KEY
    python scripts/download_entrez.py --taxon Hesionidae --email you@example.com --api-key KEY

    # By taxon ID
    python scripts/download_entrez.py --taxon 6348 --email you@example.com --api-key KEY

    # Resume interrupted download
    python scripts/download_entrez.py --taxon Phyllodocida --email you@example.com --api-key KEY --resume

Output:
    batch_001.gb, batch_002.gb, ... in output directory
    download_progress.json for resume support
"""

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

from Bio import Entrez


def parse_args():
    p = argparse.ArgumentParser(
        description="Download GenBank (.gb) files by taxon from NCBI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python scripts/download_entrez.py --taxon Phyllodocida -e you@example.com -k API_KEY
  python scripts/download_entrez.py --taxon Hesionidae -e you@example.com -k API_KEY -o hesionidae_gb
  python scripts/download_entrez.py --taxon 6348 -e you@example.com -k API_KEY --resume
  python scripts/download_entrez.py --taxon Nereis -e you@example.com -k API_KEY --validate
""",
    )
    p.add_argument("--taxon", required=True,
                   help="Taxon name (e.g. Phyllodocida, Hesionidae, Nereis) or NCBI Taxon ID (e.g. 6348)")
    p.add_argument("-o", "--output", default=None,
                   help="Output directory (default: <taxon_name>_gb/)")
    p.add_argument("-b", "--batch-size", type=int, default=500,
                   help="Records per batch file (default: 500)")
    p.add_argument("-e", "--email", required=True,
                   help="Email for NCBI Entrez (required by NCBI)")
    p.add_argument("-k", "--api-key", default=None,
                   help="NCBI API key (10 req/s, without key: 3 req/s)")
    p.add_argument("--resume", action="store_true",
                   help="Resume interrupted download from progress file")
    p.add_argument("--validate", action="store_true",
                   help="Validate downloaded files with Biopython after download")
    return p.parse_args()


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


def get_lineage(taxon_id: int) -> list[tuple[str, str]]:
    """Get taxonomy lineage as [(rank, name), ...]."""
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        lineage = []
        for entry in records[0].get("LineageEx", []):
            rank = entry.get("Rank", "no rank")
            name = entry.get("ScientificName", "")
            if rank != "no rank" and name:
                lineage.append((rank, name))
        # Add self
        rank = records[0].get("Rank", "no rank")
        name = records[0].get("ScientificName", "")
        lineage.append((rank, name))
        return lineage
    except Exception:
        return []


def resolve_taxon(taxon_input: str) -> tuple[int, str, str]:
    """Resolve taxon name or ID to (taxon_id, scientific_name, rank)."""
    if taxon_input.isdigit():
        taxon_id = int(taxon_input)
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        result = Entrez.read(handle)
        handle.close()
        info = result[0]
        return taxon_id, info["ScientificName"], info.get("Rank", "unknown")

    # Search by name
    handle = Entrez.esearch(db="taxonomy", term=taxon_input)
    result = Entrez.read(handle)
    handle.close()
    if not result["IdList"]:
        print(f"Error: taxon '{taxon_input}' not found in NCBI taxonomy", file=sys.stderr)
        sys.exit(1)

    # If multiple matches, show them
    if len(result["IdList"]) > 1:
        ids = result["IdList"][:5]
        handle = Entrez.efetch(db="taxonomy", id=",".join(ids), retmode="xml")
        summaries = Entrez.read(handle)
        handle.close()
        print(f"Multiple matches for '{taxon_input}':", file=sys.stderr)
        for s in summaries:
            print(f"  txid{s['TaxId']}: {s['ScientificName']} ({s.get('Rank', '?')})", file=sys.stderr)
        # Use first match
        taxon_id = int(result["IdList"][0])
        chosen = summaries[0]
        print(f"Using: txid{taxon_id} {chosen['ScientificName']}", file=sys.stderr)
        return taxon_id, chosen["ScientificName"], chosen.get("Rank", "unknown")

    taxon_id = int(result["IdList"][0])
    handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
    info = Entrez.read(handle)[0]
    handle.close()
    return taxon_id, info["ScientificName"], info.get("Rank", "unknown")


def search_records(taxon_id: int) -> tuple[int, str, str]:
    """Run esearch on nuccore, return (count, webenv, query_key)."""
    handle = Entrez.esearch(
        db="nuccore",
        term=f"txid{taxon_id}[Organism]",
        usehistory="y",
        retmax=0,
    )
    result = Entrez.read(handle)
    handle.close()
    count = int(result["Count"])
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]
    return count, webenv, query_key


def download_batch(webenv: str, query_key: str, retstart: int,
                   batch_size: int, output_path: Path) -> bool:
    """Download one batch of records. Returns True on success."""
    try:
        handle = Entrez.efetch(
            db="nuccore",
            rettype="gb",
            retmode="text",
            retstart=retstart,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
        )
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


def validate_files(output_dir: Path) -> tuple[int, int]:
    """Validate all .gb files. Returns (valid_count, total_records)."""
    from Bio import SeqIO
    gb_files = sorted(output_dir.glob("batch_*.gb"))
    valid = 0
    total_records = 0
    for f in gb_files:
        try:
            recs = list(SeqIO.parse(str(f), "genbank"))
            if recs:
                valid += 1
                total_records += len(recs)
            else:
                print(f"  WARNING: {f.name} has no records")
        except Exception as e:
            print(f"  INVALID: {f.name}: {e}")
    return valid, total_records


def format_eta(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        return f"{seconds / 60:.1f}min"
    return f"{seconds / 3600:.1f}h"


def main():
    args = parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # Resolve taxon
    taxon_id, taxon_name, rank = resolve_taxon(args.taxon)

    # Default output dir: <taxon_name>_gb/
    output_dir = Path(args.output) if args.output else Path(f"{taxon_name.replace(' ', '_')}_gb")
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"=" * 60)
    print(f"  Taxon: {taxon_name}")
    print(f"  Rank:  {rank}")
    print(f"  ID:    txid{taxon_id}")
    print(f"  Output: {output_dir.resolve()}")

    # Show lineage
    lineage = get_lineage(taxon_id)
    if lineage:
        lineage_str = " > ".join(name for _, name in lineage)
        print(f"  Lineage: {lineage_str}")

    # Resume or fresh search
    progress = load_progress(output_dir)
    if args.resume and progress:
        total = progress["total_count"]
        webenv = progress["webenv"]
        query_key = progress["query_key"]
        start_batch = progress.get("last_batch", 0) + 1
        print(f"  Status: Resuming from batch {start_batch}")
    else:
        total, webenv, query_key = search_records(taxon_id)
        start_batch = 0

    print(f"  Records: {total:,}")
    print(f"=" * 60)

    if total == 0:
        print("No GenBank records found. Exiting.")
        return

    num_batches = (total + args.batch_size - 1) // args.batch_size
    delay = 0.1 if args.api_key else 0.34
    eta_seconds = (num_batches - start_batch) * (delay + 0.5)  # rough estimate
    print(f"  Batches: {num_batches} ({args.batch_size}/batch)")
    print(f"  Rate limit: {'10 req/s (API key)' if args.api_key else '3 req/s (no key)'}")
    if start_batch == 0:
        print(f"  Estimated time: {format_eta(eta_seconds)}")
    print()

    # Download loop
    t_start = time.time()
    failed_count = 0

    for batch_idx in range(start_batch, num_batches):
        retstart = batch_idx * args.batch_size
        batch_file = output_dir / f"batch_{batch_idx + 1:04d}.gb"
        end = min(retstart + args.batch_size, total)

        elapsed = time.time() - t_start
        if batch_idx > start_batch:
            rate = (batch_idx - start_batch) / elapsed
            remaining = (num_batches - batch_idx) / rate if rate > 0 else 0
            eta_str = format_eta(remaining)
        else:
            eta_str = "..."
        pct = (batch_idx + 1) / num_batches * 100

        print(f"  [{pct:5.1f}%] Batch {batch_idx + 1}/{num_batches} "
              f"({retstart + 1}-{end}/{total:,}) ETA {eta_str} ... ",
              end="", flush=True)

        success = download_batch(webenv, query_key, retstart, args.batch_size, batch_file)
        if success:
            size_kb = batch_file.stat().st_size / 1024
            print(f"OK ({size_kb:.0f} KB)")
        else:
            print("FAILED")
            failed_count += 1
            if failed_count >= 3:
                print("\nToo many failures. Saving progress for --resume.")
                progress.update({
                    "taxon_id": taxon_id,
                    "taxon_name": taxon_name,
                    "rank": rank,
                    "total_count": total,
                    "webenv": webenv,
                    "query_key": query_key,
                    "last_batch": batch_idx - 1,
                    "completed": False,
                })
                save_progress(output_dir, progress)
                sys.exit(1)
            continue

        # Save progress
        progress.update({
            "taxon_id": taxon_id,
            "taxon_name": taxon_name,
            "rank": rank,
            "total_count": total,
            "webenv": webenv,
            "query_key": query_key,
            "last_batch": batch_idx,
            "completed": batch_idx == num_batches - 1,
        })
        save_progress(output_dir, progress)
        time.sleep(delay)

    # Summary
    elapsed = time.time() - t_start
    gb_files = list(output_dir.glob("batch_*.gb"))
    total_size = sum(f.stat().st_size for f in gb_files) / (1024 * 1024)
    print(f"\n{'=' * 60}")
    print(f"  Download complete!")
    print(f"  Files:  {len(gb_files)} batch files")
    print(f"  Size:   {total_size:.1f} MB")
    print(f"  Time:   {format_eta(elapsed)}")
    print(f"{'=' * 60}")

    # Optional validation
    if args.validate:
        print("\nValidating files...")
        valid, total_records = validate_files(output_dir)
        print(f"  Valid batches: {valid}/{len(gb_files)}")
        print(f"  Total records: {total_records:,}")


if __name__ == "__main__":
    main()
