"""g2t-litmethods — fill missing sequencing / assembly methods from the paper that cites the accession.

For records of a ``g2t-mitostats`` table whose platform is Unknown or whose assembler is "Not stated" /
"No assembly info in record", the accession is searched in Europe PMC (full-text search, which indexes
accession numbers quoted in the article; no Google scraping). Each hit is listed (title, DOI, PMCID) so
users can look the paper up themselves. For open-access hits the full text is downloaded, the reference
list is cut off, and the Methods text is run through the same regex rules as ``mito_methods.yaml``.
Results are *suggestions* (``lit_*`` columns), never written over the GenBank-derived columns.
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from g2t import mitostats as ms

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
_MISSING_ASM = {ms.NOT_STATED, "No assembly info in record"}


def _get(url: str, retries: int = 3) -> bytes:
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 (fixed https host)
                return r.read()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(2 * (i + 1))
    return b""


class Cache:
    """JSON file cache so reruns do not query Europe PMC again."""

    def __init__(self, path: Path):
        self.path = path
        self.d: dict[str, Any] = json.loads(path.read_text()) if path.exists() else {}

    def get(self, key: str, fn):
        if key not in self.d:
            self.d[key] = fn()
            self.path.write_text(json.dumps(self.d, ensure_ascii=False))
        return self.d[key]


def search_accession(acc: str) -> list[dict[str, str]]:
    """Europe PMC articles that mention the accession (version stripped) in their text."""
    base = acc.split(".")[0]
    q = urllib.parse.quote(f'"{base}"')
    url = f"{EPMC}/search?query={q}&format=json&pageSize=10&resultType=lite"
    res = json.loads(_get(url)).get("resultList", {}).get("result", [])
    return [{k: str(r.get(k, "")) for k in ("pmid", "pmcid", "doi", "title", "authorString", "journalTitle", "pubYear")}
            for r in res]


def methods_text(pmcid: str) -> str:
    """Plain text of an open-access article without its reference list (citations of other papers' tools
    would otherwise count as methods used)."""
    xml = _get(f"{EPMC}/{pmcid}/fullTextXML").decode("utf8", "replace")
    xml = re.sub(r"<ref-list.*?</ref-list>", " ", xml, flags=re.S)
    xml = re.sub(r"<(table-wrap|fig)\b.*?</\1>", " ", xml, flags=re.S)
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml))


def parse_methods(text: str) -> dict[str, str]:
    """Platform category/names and assemblers named in the article text (regex rules of mitostats)."""
    plats = ms.sequencing_platforms(text)
    names, _ = ms.assemblers(text)
    return {
        "lit_seq_category": ms.sequencing_category(text) if plats else "",
        "lit_seq_platforms": "; ".join(plats),
        "lit_assemblers": "; ".join(names),
    }


def needs_lookup(row: pd.Series) -> bool:
    return row.get("seq_category", "") in ("", "Unknown") or row.get("assembler_group", "") in _MISSING_ASM


def run(table: str, out: str, cache_path: str | None = None, limit: int = 0, delay: float = 0.3) -> pd.DataFrame:
    t = pd.read_csv(table, dtype=str).fillna("")
    todo = t[t.apply(needs_lookup, axis=1)]
    if limit:
        todo = todo.head(limit)
    cache = Cache(Path(cache_path or Path(out).with_suffix(".cache.json")))
    rows = []
    for _, r in todo.iterrows():
        acc = str(r["accession"])
        try:
            hits = cache.get("s:" + acc, lambda a=acc: search_accession(a))
        except Exception as e:  # network failure must not stop the batch
            rows.append({"accession": acc, "organism": r.get("organism", ""), "lit_note": f"search failed: {e}"})
            continue
        time.sleep(delay)
        if not hits:
            rows.append({"accession": acc, "organism": r.get("organism", ""), "lit_note": "no article mentions it"})
            continue
        for h in hits:
            row = {"accession": acc, "organism": r.get("organism", ""), **{f"lit_{k}": v for k, v in h.items()},
                   "lit_note": "paper listed; full text not open access"}
            if h["pmcid"]:
                try:
                    txt = cache.get("t:" + h["pmcid"], lambda p=h["pmcid"]: methods_text(p))
                    row.update(parse_methods(txt))
                    row["lit_note"] = "parsed from open-access full text; confirm by hand"
                except Exception as e:
                    row["lit_note"] = f"full text failed: {e}"
                time.sleep(delay)
            rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(out, index=False)
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="g2t-litmethods", description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("table", help="mitostats.csv from g2t-mitostats")
    ap.add_argument("-o", "--out", default="litmethods.csv")
    ap.add_argument("--limit", type=int, default=0, help="only the first N records needing a lookup (pilot)")
    ap.add_argument("--cache", help="JSON cache file (default: <out>.cache.json)")
    a = ap.parse_args(argv)
    res = run(a.table, a.out, a.cache, a.limit)
    n = res["accession"].nunique() if len(res) else 0
    print(f"{n} records looked up, {len(res)} paper rows -> {a.out}")
    if len(res) and "lit_note" in res:
        print(res["lit_note"].str.split(":").str[0].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
