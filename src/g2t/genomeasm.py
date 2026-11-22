"""Mitogenome accession -> the nuclear genome assembly of the same specimen (g2t-genome).

GenBank and ENA have no structured mitogenome-assembler field, but many mitogenomes (Darwin Tree of Life, ERGA, ...)
are by-products of a whole-genome project on the same specimen. For every record's BioSample the genome assemblies
ENA holds are looked up and their method (NCBI ``assembly_info.assembly_method`` / ``sequencing_tech``) is reported.
This describes the NUCLEAR assembly of the specimen, not the mitogenome assembler; it is a lead / context column
and is never written into the ``assembler`` column.
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from g2t.litmethods import CACHE_VERSION, ENA, Cache, _get

NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"
_COLS = ["accession", "organism", "biosample", "genome_assembly", "genome_name", "genome_level", "genome_method",
         "genome_seq_tech", "note"]


def _chunks(xs: list[str], n: int) -> list[list[str]]:
    return [xs[i:i + n] for i in range(0, len(xs), n)]


def _ena_batch(smps: list[str]) -> dict[str, list[str]]:
    q = " OR ".join(f"sample_accession={s}" for s in smps)
    url = (f"{ENA}/search?result=assembly&query={urllib.parse.quote(q)}"
           "&fields=accession,sample_accession&format=tsv&limit=0")
    out: dict[str, list[str]] = {s: [] for s in smps}
    for line in _get(url).decode().splitlines()[1:]:
        gca, _, smp = line.strip().partition("\t")
        if gca.startswith("GC") and smp in out:
            out[smp].append(gca)
    return out


def _method(m: str) -> str:
    """NCBI writes "various" when a submitter gave no method: not information."""
    return "" if m.strip().lower() == "various" else m


def _ncbi_batch(gcas: list[str]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for d in json.loads(_get(f"{NCBI}/{','.join(gcas)}/dataset_report")).get("reports", []):
        ai = d.get("assembly_info", {})
        out[d.get("accession", "").split(".")[0]] = {
            "name": ai.get("assembly_name", ""), "level": ai.get("assembly_level", ""),
            "method": _method(ai.get("assembly_method", "")), "tech": ai.get("sequencing_tech", "")}
    return out


def _fill(cache: Cache, prefix: str, keys: list[str], size: int, fn, workers: int, ttl: float | None = None) -> None:
    """Fetch the keys missing from the cache in batches, in parallel; the cache is written once per round."""
    ck = lambda k: f"{CACHE_VERSION}:{prefix}:{k}"  # noqa: E731
    ts = cache.d.setdefault("_ts", {})
    now = time.time()
    todo = [k for k in dict.fromkeys(keys) if ck(k) not in cache.d or (ttl and now - ts.get(ck(k), 0) > ttl)]
    if not todo:
        return
    batches = _chunks(todo, size)
    done = 0
    with ThreadPoolExecutor(workers) as ex:
        for b, res in zip(batches, ex.map(lambda b: _safe(fn, b), batches)):
            for k in b:
                if res is not None:  # a failed batch is not cached, so the next run retries it
                    cache.d[ck(k)] = res.get(k, [] if prefix == "ga" else {})
                    ts[ck(k)] = now
            done += len(b)
            cache.save()
            print(f"  {prefix}: {done}/{len(todo)}", flush=True)


def _safe(fn, b):
    try:
        return fn(b)
    except Exception as e:  # noqa: BLE001
        print(f"warning: batch of {len(b)} failed ({e}); will retry next run", flush=True)
        return None


def build(table: str | pd.DataFrame, cache: Cache, limit: int = 0, workers: int = 6) -> pd.DataFrame:
    t = pd.read_csv(table, dtype=str).fillna("") if not isinstance(table, pd.DataFrame) else table.fillna("")
    t = t[t.get("biosample", "") != ""]
    if limit:
        t = t.head(limit)
    smps = [s for b in t["biosample"] for s in b.split("; ") if s]
    _fill(cache, "ga", smps, 25, _ena_batch, workers, ttl=30 * 86400)
    gcas = sorted({g.split(".")[0] for s in smps for g in cache.d.get(f"{CACHE_VERSION}:ga:{s}", [])})
    _fill(cache, "gi", gcas, 20, _ncbi_batch, workers)
    rows: list[dict[str, str]] = []
    for _, r in t.iterrows():
        gl = [g for s in r["biosample"].split("; ") if s for g in cache.d.get(f"{CACHE_VERSION}:ga:{s}", [])]
        base = {"accession": r["accession"], "organism": r.get("organism", ""), "biosample": r["biosample"]}
        if not gl:
            rows.append({**base, "note": "no genome assembly for this BioSample in ENA"})
        for g in gl:
            d = cache.d.get(f"{CACHE_VERSION}:gi:{g.split('.')[0]}", {})
            rows.append({**base, "genome_assembly": g, "genome_name": d.get("name", ""),
                         "genome_level": d.get("level", ""), "genome_method": d.get("method", ""),
                         "genome_seq_tech": d.get("tech", ""), "note": "" if d else "no NCBI assembly report"})
    return pd.DataFrame(rows, columns=_COLS)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="g2t-genome", description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("table", help="mitostats.csv from g2t-mitostats")
    ap.add_argument("-o", "--out", default="genome_assemblies.csv")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6, help="parallel requests (default 6)")
    ap.add_argument("--cache", help="JSON cache (default: <out>.cache.json)")
    a = ap.parse_args(argv)
    res = build(a.table, Cache(Path(a.cache or Path(a.out).with_suffix(".cache.json"))), a.limit, a.workers)
    res.to_csv(a.out, index=False)
    ok = res[res["genome_assembly"] != ""]
    print(f"{res.accession.nunique()} records, {ok.accession.nunique()} with a genome assembly, "
          f"{(ok.genome_method != '').sum()} rows with method -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
