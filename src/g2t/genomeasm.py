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
import urllib.parse
from pathlib import Path

import pandas as pd

from g2t.litmethods import CACHE_VERSION, ENA, Cache, _get

NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"
_COLS = ["accession", "organism", "biosample", "genome_assembly", "genome_name", "genome_level", "genome_method",
         "genome_seq_tech", "note"]


def biosample_assemblies(smp: str, cache: Cache) -> list[str]:
    def q() -> list[str]:
        url = (f"{ENA}/search?result=assembly&query={urllib.parse.quote(f'sample_accession={smp}')}"
               "&fields=accession&format=tsv&limit=0")
        return [x.strip() for x in _get(url).decode().splitlines()[1:] if x.strip().startswith("GC")]
    return cache.get(f"ga:{smp}", q, ttl=30 * 86400)


def assembly_info(gcas: list[str], cache: Cache) -> dict[str, dict[str, str]]:
    gcas = [g.split('.')[0] for g in gcas]
    todo = [g for g in gcas if f"{CACHE_VERSION}:gi:{g}" not in cache.d]
    for i in range(0, len(todo), 20):
        part = todo[i:i + 20]
        rep = json.loads(_get(f"{NCBI}/{','.join(part)}/dataset_report")).get("reports", [])
        for d in rep:
            ai = d.get("assembly_info", {})
            val = {"name": ai.get("assembly_name", ""), "level": ai.get("assembly_level", ""),
                   "method": ai.get("assembly_method", ""), "tech": ai.get("sequencing_tech", "")}
            cache.get(f"gi:{d.get('accession', '').split('.')[0]}", lambda v=val: v)
    return {g: cache.d[f"{CACHE_VERSION}:gi:{g}"] for g in gcas if f"{CACHE_VERSION}:gi:{g}" in cache.d}


def build(table: str | pd.DataFrame, cache: Cache, limit: int = 0) -> pd.DataFrame:
    t = pd.read_csv(table, dtype=str).fillna("") if not isinstance(table, pd.DataFrame) else table.fillna("")
    t = t[t.get("biosample", "") != ""]
    if limit:
        t = t.head(limit)
    rows: list[dict[str, str]] = []
    cache.defer = True
    smp_cache: dict[str, list[str]] = {}
    for n, (_, r) in enumerate(t.iterrows(), 1):
        if n % 100 == 0:
            cache.save()
        gcas: list[str] = []
        for smp in [x for x in r["biosample"].split("; ") if x]:
            smp_cache.setdefault(smp, biosample_assemblies(smp, cache))
            gcas += smp_cache[smp]
        info = assembly_info(gcas, cache) if gcas else {}
        base = {"accession": r["accession"], "organism": r.get("organism", ""), "biosample": r["biosample"]}
        if not gcas:
            rows.append({**base, "note": "no genome assembly for this BioSample in ENA"})
        for g in gcas:
            d = info.get(g, {})
            rows.append({**base, "genome_assembly": g, "genome_name": d.get("name", ""),
                         "genome_level": d.get("level", ""), "genome_method": d.get("method", ""),
                         "genome_seq_tech": d.get("tech", ""),
                         "note": "" if d else "no NCBI assembly report"})
    cache.save()
    cache.defer = False
    return pd.DataFrame(rows, columns=_COLS)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="g2t-genome", description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("table", help="mitostats.csv from g2t-mitostats")
    ap.add_argument("-o", "--out", default="genome_assemblies.csv")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--cache", help="JSON cache (default: <out>.cache.json)")
    a = ap.parse_args(argv)
    res = build(a.table, Cache(Path(a.cache or Path(a.out).with_suffix(".cache.json"))), a.limit)
    res.to_csv(a.out, index=False)
    ok = res[res["genome_assembly"] != ""]
    print(f"{res.accession.nunique()} records, {ok.accession.nunique()} with a genome assembly, "
          f"{(ok.genome_method != '').sum()} rows with method -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
