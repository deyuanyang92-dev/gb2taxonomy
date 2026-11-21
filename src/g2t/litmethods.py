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
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pandas as pd

from g2t import mitostats as ms

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
_MISSING_ASM = {ms.NOT_STATED, "No assembly info in record"}
CACHE_VERSION = "v2"  # bump when extraction rules change; raw XML stays cached, parsing reruns every time
MAX_HITS = 100
ENA = "https://www.ebi.ac.uk/ena/portal/api"


def _get(url: str, retries: int = 4) -> bytes:
    """GET with retries for network errors, HTTP 429/5xx (honours Retry-After); other HTTP errors are permanent."""
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 (fixed https host)
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or i == retries - 1:
                raise
            wait = e.headers.get("Retry-After", "") if e.headers else ""
            time.sleep(min(int(wait), 120) if wait.isdigit() else 2 * (i + 1))
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(2 * (i + 1))
    return b""


SEARCH_TTL = 30 * 86400  # Europe PMC search hits are re-queried after 30 days


class Cache:
    """JSON file cache so reruns do not query Europe PMC again; written atomically, a damaged file is ignored."""

    def __init__(self, path: Path):
        self.path = path
        self.d: dict[str, Any] = {}
        if path.exists():
            try:
                self.d = json.loads(path.read_text())
            except ValueError:
                print(f"warning: cache {path} is damaged, starting empty")

    def get(self, key: str, fn, ttl: float | None = None):
        """Cached fn(); with ttl (seconds) an entry older than that is fetched again (search results go stale:
        a paper may be indexed later)."""
        key = f"{CACHE_VERSION}:{key}" if not key.startswith("x:") else key
        ts = self.d.setdefault("_ts", {})
        if key not in self.d or (ttl is not None and time.time() - ts.get(key, 0) > ttl):
            self.d[key] = fn()
            ts[key] = time.time()
            tmp = self.path.with_name(self.path.name + ".tmp")
            tmp.write_text(json.dumps(self.d, ensure_ascii=False))
            os.replace(tmp, self.path)
        return self.d[key]


def search_accession(acc: str) -> list[dict[str, str]]:
    """Europe PMC articles that mention the accession (version stripped) in their text (cursor-paged, <= MAX_HITS)."""
    base = acc.split(".")[0]
    q = urllib.parse.quote(f'"{base}"')
    out: list[dict[str, str]] = []
    cursor = "*"
    while len(out) < MAX_HITS:
        url = f"{EPMC}/search?query={q}&format=json&pageSize=25&resultType=lite&cursorMark={urllib.parse.quote(cursor)}"
        data = json.loads(_get(url))
        res = data.get("resultList", {}).get("result", [])
        out += [{k: str(r.get(k, "")) for k in ("pmid", "pmcid", "doi", "title", "authorString", "journalTitle",
                                                "pubYear")} for r in res]
        nxt = data.get("nextCursorMark", "")
        if not res or not nxt or nxt == cursor:
            break
        cursor = nxt
    return out


def search_pmid(pmid: str) -> list[dict[str, str]]:
    """The article a record's own REFERENCE ... PUBMED field points to (Europe PMC, resultType lite)."""
    q = urllib.parse.quote(f"EXT_ID:{pmid} AND SRC:MED")
    data = json.loads(_get(f"{EPMC}/search?query={q}&format=json&pageSize=5&resultType=lite"))
    return [{k: str(r.get(k, "")) for k in ("pmid", "pmcid", "doi", "title", "authorString", "journalTitle",
                                            "pubYear")} for r in data.get("resultList", {}).get("result", [])]


def fetch_xml(pmcid: str) -> str:
    return _get(f"{EPMC}/{pmcid}/fullTextXML").decode("utf8", "replace")


_METHOD_TITLE = re.compile(r"method|material|procedure|sequenc|assembl|experimental|sampling", re.I)
_TEXT = re.compile(r"\s+")


def _text(el: ET.Element) -> str:
    return _TEXT.sub(" ", " ".join(el.itertext())).strip()


def article_sections(xml: str) -> tuple[str, str, bool]:
    """(Methods text, whole-body text, methods_found). Reference list, tables and figures are dropped; Methods are
    the body sections whose sec-type or title says so (nested subsections included)."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        txt = _TEXT.sub(" ", re.sub(r"<[^>]+>", " ", re.sub(r"<ref-list.*?</ref-list>", " ", xml, flags=re.S)))
        return txt, txt, False
    for parent in root.iter():
        for ch in list(parent):
            if ch.tag in ("ref-list", "table-wrap", "fig", "back"):
                parent.remove(ch)
    body = root.find(".//body")
    if body is None:
        body = root
    meth = []
    for sec in body.iter("sec"):
        title = sec.findtext("title") or ""
        if _METHOD_TITLE.search(title) or _METHOD_TITLE.search(sec.get("sec-type", "")):
            meth.append(_text(sec))
    whole = _text(body)
    return (" ".join(meth), whole, True) if meth else (whole, whole, False)


def role_text(xml: str) -> str:
    """Whole article text for the reports-vs-cited decision: only the reference list is dropped, so deposits stated in
    back matter (data availability), tables or figure captions are still seen."""
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return article_sections(xml)[1]
    for parent in root.iter():
        for ch in list(parent):
            if ch.tag == "ref-list":
                parent.remove(ch)
    return _text(root)


def accession_re(acc: str) -> re.Pattern:
    """Exact accession token: AB123456 matches AB123456 / AB123456.1, never AB1234567 or XAB123456."""
    return re.compile(r"(?<![A-Za-z0-9])" + re.escape(acc.split(".")[0]) + r"(?:\.\d+)?(?![A-Za-z0-9])")


_SENT = re.compile(r"(?<=[.;])\s+(?=[A-Z(])")
_DEPOSIT = re.compile(r"submitted|deposit|accession (?:no|number)|registered|available (?:in|at|from|under)|GenBank under|"
                      r"under (?:the )?accession|has been (?:uploaded|released)|were? (?:uploaded|released)", re.I)
_REUSE = re.compile(r"download|retriev|obtain\w* from|taken from|sequences? (?:of|from) (?:other|related)|"
                    r"outgroup|used in (?:the )?(?:phylogen|analys)|reference sequence", re.I)


def accession_role(body: str, acc: str) -> tuple[str, str]:
    """("reports" | "cited", evidence sentence): does the article say it *deposited* the sequence, or only mention it?"""
    rx = accession_re(acc)
    # a character window, not a sentence: "accession no. OP104125" would be cut at "no."
    wins = [body[max(0, m.start() - 220): m.end() + 60] for m in rx.finditer(body)]
    for w in wins:
        if _DEPOSIT.search(w) and not _REUSE.search(w):
            return "reports", w
    return "cited", (wins[0] if wins else "")


_SEQ_CTX = re.compile(r"sequenc|librar|reads?\b|platform|paired[- ]end|flow ?cell|nanopore|hifi|run\b", re.I)
_ASM_CTX = re.compile(r"assembl|de novo|reconstruct|circulari", re.I)
_ANNOT = re.compile(r"\bannotat\w*", re.I)
_NOT_ASM = re.compile(r"alignment|aligned|phylogen|nuclear|primer", re.I)


def parse_methods(text: str, mito_paper: bool = True) -> dict[str, str]:
    """Platform / assemblers named in sequencing- and assembly-context sentences (not in the whole text), with evidence."""
    sents = [x for x in _SENT.split(text) if x]
    plat_s = [x for x in sents if _SEQ_CTX.search(x)]
    plats = ms.sequencing_platforms(" ".join(plat_s))
    names: list[str] = []
    ev = []
    for x0 in sents:
        # "assembled with MitoZ and annotated with MITOS": only the part before "annotat…" can name the assembler
        x = _ANNOT.split(x0)[0]
        if _ASM_CTX.search(x) and not _NOT_ASM.search(x) and (mito_paper or ms._MITO_WORD.search(x)):
            hit = [n for n, rx in ms._rules()["assemblers"] if rx.search(x) and n not in names]
            if hit:
                names += hit
                ev.append(x[:200])
    return {
        "lit_seq_category": ms.sequencing_category(" ".join(plat_s)) if plats else "",
        "lit_seq_platforms": "; ".join(plats),
        "lit_assemblers": "; ".join(names),
        "lit_evidence": " | ".join(ev),
    }


_GENOMIC = {"wgs", "wga", "wxs"}
_PLATFORM_WORD = {"PACBIO_SMRT": "PacBio", "OXFORD_NANOPORE": "Nanopore", "LS454": "454", "ION_TORRENT": "Ion Torrent"}


def ena_runs(runs: list[str], cache: Cache) -> dict[str, dict[str, str]]:
    """ENA read_run metadata (instrument platform / model, library strategy) for run accessions; INSDC mirrors
    NCBI SRA runs (SRR/DRR) too. Queried 50 at a time and cached per run."""
    todo = [r for r in runs if f"{CACHE_VERSION}:r:{r}" not in cache.d]
    for i in range(0, len(todo), 50):
        part = todo[i:i + 50]
        q = urllib.parse.quote(" OR ".join(f"run_accession={r}" for r in part))
        url = (f"{ENA}/search?result=read_run&query={q}&fields=run_accession,instrument_platform,instrument_model,"
               "library_strategy&format=tsv&limit=0")
        got: dict[str, dict[str, str]] = {}
        for line in _get(url).decode().splitlines()[1:]:
            f = line.split("\t") + [""] * 4
            got[f[0]] = {"platform": f[1], "model": f[2], "strategy": f[3]}
        for r in part:  # only runs ENA returned are cached: an unknown run may appear later and is asked again
            if r in got:
                cache.get(f"r:{r}", lambda r=r: got[r])
    return {r: cache.d[f"{CACHE_VERSION}:r:{r}"] for r in runs if f"{CACHE_VERSION}:r:{r}" in cache.d}


def ena_sample_runs(samples: list[str], cache: Cache) -> list[str]:
    """Run accessions ENA holds for BioSamples (the specimen the mitogenome came from; a study / BioProject would mix
    specimens and is not used). Found runs are cached with their instrument fields."""
    runs: list[str] = []
    for smp in samples:
        def fetch(smp=smp):
            q = urllib.parse.quote(f"sample_accession={smp} OR secondary_sample_accession={smp}")
            url = (f"{ENA}/search?result=read_run&query={q}&fields=run_accession,instrument_platform,instrument_model,"
                   "library_strategy&format=tsv&limit=0")
            lines = _get(url).decode().splitlines()[1:]
            return [ln.split("\t") + [""] * 4 for ln in lines if ln.strip()]
        for f in cache.get(f"b:{smp}", fetch, ttl=SEARCH_TTL):
            cache.get(f"r:{f[0]}", lambda f=f: {"platform": f[1], "model": f[2], "strategy": f[3]})
            runs.append(f[0])
    return list(dict.fromkeys(runs))


def sra_methods(runs: list[str], cache: Cache) -> dict[str, str]:
    """Sequencing platform of a record from its raw-read runs (instrument fields, no text matching). Only runs whose
    library strategy is positively genomic (WGS / WGA / WXS) count; Hi-C, RNA-seq and ambiguous ones (OTHER, empty)
    are ignored, and with no genomic run nothing is inferred. Runs ENA does not return are listed in
    `sra_unresolved`: the platform is then only a partial answer (`sra_complete` = no)."""
    meta = ena_runs(runs, cache)
    unresolved = [r for r in runs if r not in meta]
    genomic = {r: m for r, m in meta.items() if m["strategy"].lower() in _GENOMIC}
    if not genomic:
        return {"sra_unresolved": "; ".join(unresolved)} if unresolved else {}
    text = " ".join(f"{_PLATFORM_WORD.get(m['platform'], m['platform'])} {m['model']}" for m in genomic.values())
    plats = ms.sequencing_platforms(text)
    return {
        "sra_seq_category": ms.sequencing_category(text) if plats else "",
        "sra_seq_platforms": "; ".join(plats),
        "sra_instruments": "; ".join(sorted({m["model"] for m in genomic.values() if m["model"]})),
        "sra_strategies": "; ".join(sorted({m["strategy"] for m in meta.values() if m["strategy"]})),
        "sra_runs_used": "; ".join(genomic),
        "sra_unresolved": "; ".join(unresolved),
        "sra_complete": "no" if unresolved else "yes",
    }


def needs_lookup(row: pd.Series) -> bool:
    # a platform already named in the record (COMMENT, Assembly-Data) is not overridden: only unknown or the bare "NGS"
    # category (SRA link without instrument) is looked up in ENA
    no_platform = row.get("seq_category", "") in ("", "Unknown") or (
        (row.get("sra_runs", "") or row.get("biosample", "")) and row.get("seq_category", "") == "NGS" and not row.get("seq_platforms", ""))
    return bool(no_platform) or row.get("assembler_group", "") in _MISSING_ASM


def run(table: str, out: str, cache_path: str | None = None, limit: int = 0, delay: float = 0.3) -> pd.DataFrame:
    t = pd.read_csv(table, dtype=str).fillna("")
    todo = t[t.apply(needs_lookup, axis=1)]
    if limit:
        todo = todo.head(limit)
    cache = Cache(Path(cache_path or Path(out).with_suffix(".cache.json")))
    rows = []
    for _, r in todo.iterrows():
        acc = str(r["accession"])
        base = {"accession": acc, "organism": r.get("organism", "")}
        if (r.get("sra_runs") or r.get("biosample")) and r.get("seq_category", "") in ("", "Unknown", "NGS") \
                and not r.get("seq_platforms", ""):
            try:
                runs = r["sra_runs"].split("; ") if r.get("sra_runs") else []
                via = "ENA read_run instrument fields"
                if not runs:
                    runs = ena_sample_runs(r["biosample"].split("; "), cache)
                    via = "ENA read_run of the record's BioSample"
                if runs:
                    base.update(sra_methods(runs, cache))
                if base.get("sra_seq_platforms"):
                    base["method_source"] = via + ("" if base.get("sra_complete") == "yes" else " (partial run set)")
            except Exception as e:
                base["sra_note"] = f"ENA read_run failed: {e}"
        try:
            hits = cache.get("s:" + acc, lambda a=acc: search_accession(a), ttl=SEARCH_TTL)
        except Exception as e:  # network failure must not stop the batch
            rows.append({**base, "lit_note": f"search failed: {e}"})
            continue
        time.sleep(delay)
        own: set[str] = set()  # papers the record itself cites in its REFERENCE block
        for pm in [x for x in r.get("pubmed_ids", "").split("; ") if x]:
            try:
                for h in cache.get("p:" + pm, lambda pm=pm: search_pmid(pm), ttl=SEARCH_TTL):
                    own.add(h["pmid"])
                    if all(h["pmid"] != x["pmid"] for x in hits):
                        hits = hits + [h]
            except Exception as e:
                base["lit_note_pmid"] = f"PubMed id lookup failed: {e}"
        if not hits:
            rows.append({**base, "lit_note": "no article mentions it"})
            continue
        for h in hits:
            row = {**base, **{f"lit_{k}": v for k, v in h.items()},
                   "lit_via": "record REFERENCE PUBMED" if h["pmid"] in own else "accession in article text",
                   "lit_links": "; ".join(x for x in (
                       f"https://doi.org/{h['doi']}" if h["doi"] else "",
                       f"https://pubmed.ncbi.nlm.nih.gov/{h['pmid']}/" if h["pmid"] else "",
                       f"https://europepmc.org/article/MED/{h['pmid']}" if h["pmid"] else "",
                       ) if x),
                   "lit_note": "paper listed; full text not open access"}
            if h["pmcid"]:
                try:
                    xml = cache.get("x:" + h["pmcid"], lambda p=h["pmcid"]: fetch_xml(p))
                    meth, whole, found = article_sections(xml)
                    role, ev = accession_role(role_text(xml), acc)
                    if role == "cited" and not ev and h["pmid"] in own:
                        role, ev = "record_reference", "paper is the record's own REFERENCE (PUBMED field)"
                    row["lit_role"], row["lit_accession_sentence"] = role, ev
                    if role in ("reports", "record_reference"):
                        row.update(parse_methods(meth, mito_paper=bool(re.search("mito", h["title"], re.I))))
                        row["method_source"] = row.get("method_source") or "literature Methods (Europe PMC)"
                        row["lit_note"] = ("parsed from Methods section" if found else
                                           "no Methods section found, whole text parsed") + "; confirm by hand"
                    else:
                        row["lit_note"] = "article mentions the accession but does not report depositing it; methods not taken"
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
