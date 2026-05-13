#!/usr/bin/env python3
"""
Download GenBank records by taxon using Python requests (direct NCBI E-utilities API).

Usage:
    python scripts/download_requests.py --taxon 6270 --email you@example.com --api-key YOUR_KEY
    python scripts/download_requests.py --taxon 6270 --email you@example.com --resume

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

import requests


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def parse_args():
    p = argparse.ArgumentParser(
        description="Download GenBank records by taxon (requests + E-utilities method)"
    )
    p.add_argument("--taxon", required=True, help="NCBI Taxon ID (e.g. 6270) or name")
    p.add_argument("--output", default="./gb_download", help="Output directory")
    p.add_argument("--batch-size", type=int, default=500, help="Records per batch (default: 500)")
    p.add_argument("--email", required=True, help="Email for NCBI")
    p.add_argument("--api-key", default=None, help="NCBI API key (10 req/s)")
    p.add_argument("--max-retries", type=int, default=3, help="Max retries per request")
    p.add_argument("--timeout", type=int, default=120, help="HTTP timeout in seconds")
    p.add_argument("--resume", action="store_true", help="Resume from previous download")
    return p.parse_args()


def fetch_with_retry(session: requests.Session, url: str, params: dict,
                     max_retries: int = 3, timeout: int = 120) -> requests.Response:
    """HTTP GET with exponential backoff retry."""
    for attempt in range(max_retries):
        try:
            resp = session.get(url, params=params, timeout=timeout)
            if resp.status_code == 200:
                return resp
            if resp.status_code == 429:
                wait = 10 * (2 ** attempt)
                print(f"    Rate limited, waiting {wait}s...", flush=True)
                time.sleep(wait)
                continue
            if resp.status_code >= 500:
                wait = 2 * (2 ** attempt)
                print(f"    Server error {resp.status_code}, retry in {wait}s...", flush=True)
                time.sleep(wait)
                continue
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                wait = 2 * (2 ** attempt)
                print(f"    Network error: {e}, retry in {wait}s...", flush=True)
                time.sleep(wait)
            else:
                raise
    raise RuntimeError(f"Failed after {max_retries} retries: {url}")


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


def resolve_taxon(session: requests.Session, taxon_input: str, api_key: str = None) -> tuple[int, str]:
    """Resolve taxon name or ID."""
    params = {"db": "taxonomy", "retmode": "json"}
    if api_key:
        params["api_key"] = api_key

    if taxon_input.isdigit():
        taxon_id = int(taxon_input)
        params["id"] = taxon_id
    else:
        params["term"] = taxon_input
        resp = fetch_with_retry(session, ESEARCH_URL, params)
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            print(f"Error: taxon '{taxon_input}' not found", file=sys.stderr)
            sys.exit(1)
        taxon_id = int(ids[0])
        params = {"db": "taxonomy", "id": taxon_id, "retmode": "json"}
        if api_key:
            params["api_key"] = api_key

    resp = fetch_with_retry(session, ESUMMARY_URL, params)
    data = resp.json()
    result = data.get("result", {})
    name = result.get(str(taxon_id), {}).get("scientificname", f"txid{taxon_id}")
    return taxon_id, name


def search_records(session: requests.Session, taxon_id: int,
                   api_key: str = None) -> tuple[int, str, str]:
    """esearch to get count + WebEnv + QueryKey."""
    params = {
        "db": "nuccore",
        "term": f"txid{taxon_id}[Organism]",
        "usehistory": "y",
        "retmax": 0,
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    resp = fetch_with_retry(session, ESEARCH_URL, params)
    data = resp.json()["esearchresult"]
    count = int(data["count"])
    webenv = data["webenv"]
    query_key = data["querykey"]
    return count, webenv, query_key


def download_batch(session: requests.Session, webenv: str, query_key: str,
                   retstart: int, batch_size: int, output_path: Path,
                   api_key: str = None, timeout: int = 120) -> bool:
    """Download one batch. Returns True on success."""
    params = {
        "db": "nuccore",
        "rettype": "gb",
        "retmode": "text",
        "retstart": retstart,
        "retmax": batch_size,
        "webenv": webenv,
        "query_key": query_key,
    }
    if api_key:
        params["api_key"] = api_key

    try:
        resp = fetch_with_retry(session, EFETCH_URL, params, timeout=timeout)
        text = resp.text
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

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "g2t-download/1.0"})

    # Resume: load existing progress
    progress = load_progress(output_dir)
    if args.resume and progress:
        taxon_id = progress["taxon_id"]
        taxon_name = progress["taxon_name"]
        total = progress["total_count"]
        webenv = progress["webenv"]
        query_key = progress["query_key"]
        start_batch = progress.get("last_batch", 0) + 1
        print(f"Resuming: {taxon_name} (txid{taxon_id}), {total} records, starting batch {start_batch}")
    else:
        taxon_id, taxon_name = resolve_taxon(session, args.taxon, args.api_key)
        total, webenv, query_key = search_records(session, taxon_id, args.api_key)
        start_batch = 0
        print(f"Taxon: {taxon_name} (txid{taxon_id}), {total} records found")

    if total == 0:
        print("No records found. Exiting.")
        return

    num_batches = (total + args.batch_size - 1) // args.batch_size
    print(f"Will download in {num_batches} batches of {args.batch_size} records")

    delay = 0.1 if args.api_key else 0.34

    for batch_idx in range(start_batch, num_batches):
        retstart = batch_idx * args.batch_size
        batch_file = output_dir / f"batch_{batch_idx + 1:04d}.gb"
        end = min(retstart + args.batch_size, total)

        print(f"  Batch {batch_idx + 1}/{num_batches} ({retstart + 1}-{end}/{total})...", end=" ", flush=True)

        success = download_batch(
            session, webenv, query_key, retstart, args.batch_size,
            batch_file, args.api_key, args.timeout
        )

        if success:
            size_kb = batch_file.stat().st_size / 1024
            print(f"OK ({size_kb:.1f} KB)")
        else:
            print("FAILED (will retry on --resume)")
            progress.update({
                "taxon_id": taxon_id,
                "taxon_name": taxon_name,
                "total_count": total,
                "webenv": webenv,
                "query_key": query_key,
                "last_batch": batch_idx - 1,
                "completed": False,
            })
            save_progress(output_dir, progress)
            sys.exit(1)

        progress.update({
            "taxon_id": taxon_id,
            "taxon_name": taxon_name,
            "total_count": total,
            "webenv": webenv,
            "query_key": query_key,
            "last_batch": batch_idx,
            "completed": batch_idx == num_batches - 1,
        })
        save_progress(output_dir, progress)

        time.sleep(delay)

    gb_files = list(output_dir.glob("batch_*.gb"))
    total_size = sum(f.stat().st_size for f in gb_files) / (1024 * 1024)
    print(f"\nDone! {len(gb_files)} batch files, {total_size:.1f} MB total")


if __name__ == "__main__":
    main()
