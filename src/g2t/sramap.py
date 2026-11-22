"""Mitogenome accession -> the raw SRA/ENA runs behind it (g2t-sra).

For every record of a ``g2t-mitostats`` table the runs are taken from the record itself (DBLINK / COMMENT ``sra_runs``);
a record without runs falls back to the reads ENA holds for its BioSample (the specimen; a BioProject would mix
specimens and is not used). One row per accession x run with ENA's read_run metadata, so the platform, library
strategy and read volume of the data an assembly was made from can be checked. ``run_link`` says where the run came
from (``record`` or ``biosample``). Results are lookups, never written over GenBank-derived columns.
"""
from __future__ import annotations

import argparse
import urllib.parse
from pathlib import Path

import pandas as pd

from g2t.litmethods import CACHE_VERSION, ENA, Cache, _get

FIELDS = ["run_accession", "experiment_accession", "study_accession", "sample_accession", "secondary_sample_accession",
          "instrument_platform", "instrument_model", "library_strategy", "library_source", "library_layout",
          "read_count", "base_count", "first_public"]
_NEW = ["run", "experiment", "study", "sample", "biosample", "platform", "model", "strategy", "source", "layout",
        "reads", "bases", "first_public"]
GENOMIC = {"wgs", "wga", "wxs"}


def _query(q: str) -> list[dict[str, str]]:
    url = f"{ENA}/search?result=read_run&query={urllib.parse.quote(q)}&fields={','.join(FIELDS)}&format=tsv&limit=0"
    lines = _get(url).decode().splitlines()
    return [dict(zip(_NEW, ln.split("\t") + [""] * len(_NEW))) for ln in lines[1:] if ln.strip()]


def runs_meta(runs: list[str], cache: Cache) -> dict[str, dict[str, str]]:
    todo = [r for r in runs if f"{CACHE_VERSION}:sr:{r}" not in cache.d]
    for i in range(0, len(todo), 50):
        part = todo[i:i + 50]
        got = {d["run"]: d for d in _query(" OR ".join(f"run_accession={r}" for r in part))}
        for r in part:
            if r in got:  # unknown runs are not cached: ENA may get them later
                cache.get(f"sr:{r}", lambda d=got[r]: d)
    return {r: cache.d[f"{CACHE_VERSION}:sr:{r}"] for r in runs if f"{CACHE_VERSION}:sr:{r}" in cache.d}


def biosample_runs(smp: str, cache: Cache) -> list[dict[str, str]]:
    return cache.get(f"sb:{smp}", lambda: _query(f"sample_accession={smp} OR secondary_sample_accession={smp}"),
                     ttl=30 * 86400)


def build(table: str | pd.DataFrame, cache: Cache, limit: int = 0) -> pd.DataFrame:
    t = pd.read_csv(table, dtype=str).fillna("") if not isinstance(table, pd.DataFrame) else table.fillna("")
    cols = [c for c in ("sra_runs", "biosample") if c in t]
    t = t[t[cols].apply(lambda r: any(r[c] for c in cols), axis=1)] if cols else t.iloc[0:0]
    if limit:
        t = t.head(limit)
    rows: list[dict[str, str]] = []
    cache.defer = True
    for n, (_, r) in enumerate(t.iterrows(), 1):
        if n % 100 == 0:
            cache.save()
        base = {"accession": r["accession"], "organism": r.get("organism", ""),
                "biosample_in_record": r.get("biosample", "")}
        runs = [x for x in r.get("sra_runs", "").split("; ") if x]
        link, found = "record", []
        if runs:
            meta = runs_meta(runs, cache)
            found = [meta.get(x, {"run": x}) for x in runs]
        else:
            link = "biosample"
            for smp in [x for x in r.get("biosample", "").split("; ") if x]:
                found += biosample_runs(smp, cache)
        if not found:
            rows.append({**base, "run_link": link, "note": "no run found in ENA"})
        for d in found:
            rows.append({**base, "run_link": link, **d,
                         "note": "" if d.get("strategy") else "run not returned by ENA",
                         "genomic": "yes" if d.get("strategy", "").lower() in GENOMIC else "no"})
    cache.save()
    cache.defer = False
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="g2t-sra", description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("table", help="mitostats.csv from g2t-mitostats")
    ap.add_argument("-o", "--out", default="sra_runs.csv")
    ap.add_argument("--limit", type=int, default=0, help="only the first N records with runs / BioSample (pilot)")
    ap.add_argument("--cache", help="JSON cache (default: <out>.cache.json)")
    a = ap.parse_args(argv)
    res = build(a.table, Cache(Path(a.cache or Path(a.out).with_suffix(".cache.json"))), a.limit)
    res.to_csv(a.out, index=False)
    n = res["accession"].nunique() if len(res) else 0
    ok = res[res.get("run", pd.Series(dtype=str)).fillna("") != ""] if len(res) else res
    print(f"{n} records, {len(ok)} runs ({ok['accession'].nunique() if len(ok) else 0} records with runs) -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
