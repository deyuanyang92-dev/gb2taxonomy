#!/usr/bin/env python3
"""Download the GenBank records of a taxon from NCBI nuccore (g2t step 0).

Default selection ("markers"): no WGS contigs, mRNA, RefSeq, or nuclear records > 100 kb.
Records are kept by accession.version in a SQLite store (or only in the batch files with
--no-store); a run fetches only records it does not have, --since limits the NCBI query to
records modified since the last run. Output: <out>/<tag>/batch_NNNN.gb, accessions.tsv,
manifest.json, changes.tsv.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import socket
import sys
import time
import urllib.parse
import urllib.request
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from datetime import date, timedelta
from pathlib import Path
from typing import Callable

from g2t.ncbi_genes import MITO_GENES, NUCLEAR_GENES, canon, clause
from g2t.recstore import DirStore, RecordStore, iter_records

# Nuclear records above LARGE_BP (chromosomes, scaffolds, TSA masters) are skipped unless --include-large.
# Organelle records are never length-capped. plastid[filter] matches nothing and mitochondri* is truncated by Entrez.
LARGE_BP = 100_000
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TIMEOUT = 60          # s without data before a request is retried
ORGANELLE = ("(mitochondrion[filter] OR chloroplast[filter] OR mitochondrion[Title] OR mitochondrial[Title] "
             "OR chloroplast[Title] OR plastid[Title])")
# also matches "<species> genome assembly, organelle: mitochondrion" (Darwin Tree of Life); no upper length
MITOGENOME = ('AND (mitochondrion[filter] OR mitochondrion[Title] OR mitochondrial[Title]) AND (complete genome[Title] '
              'OR "complete mitochondrial genome"[Title] OR mitogenome[Title] OR "mitochondrial genome"[Title] '
              'OR "genome assembly"[Title]) AND 10000:999999999[SLEN]')

warnings.filterwarnings("ignore", category=UserWarning, module="Bio.Entrez")


@dataclass
class DownloadOptions:
    all: bool = False
    mito: bool = False
    mitogenome: bool = False
    gene: str = ""
    minlen: int | None = None
    maxlen: int | None = None
    query: str = ""
    include_wgs: bool = False
    include_mrna: bool = False
    include_refseq: bool = False
    include_large: bool = False
    since: str = ""           # "" = full list; "auto" = last download date; or YYYY-MM-DD ([MDAT] incremental)
    no_store: bool = False    # keep records only in the selection's batch files (no SQLite store)
    batch: int = 500
    workers: int = 3
    tag: str = ""


@dataclass
class DownloadResult:
    taxid: str
    name: str
    rank: str
    tag: str
    query: str
    count: int
    out_dir: str
    composition: dict[str, int] = field(default_factory=dict)
    per_gene: dict[str, int] = field(default_factory=dict)
    batches: int = 0          # batches fetched from NCBI in this run
    fetched: int = 0          # records fetched from NCBI in this run
    reused: int = 0           # records taken from the local store
    to_fetch: int = 0         # records not in the store (dry run: would be fetched)
    failed: int = 0           # records that failed after all retries and batch splitting
    new: int = 0              # changes against the previous accession list of this selection
    updated: int = 0
    removed: int = 0
    query_changed: bool = False   # the selection's query differs from the previous run (g2t default changed)
    since: str = ""               # Entrez date of an incremental run (records modified since then)
    candidates: int = 0           # records modified since that date


def setup_entrez(email: str | None = None, api_key: str | None = None) -> float:
    """Configure Bio.Entrez; returns the delay between requests (s)."""
    from Bio import Entrez
    Entrez.email = email or os.environ.get("NCBI_EMAIL") or None  # type: ignore[assignment]
    Entrez.api_key = api_key or os.environ.get("NCBI_API_KEY") or None  # type: ignore[assignment]
    Entrez.tool = "g2t"
    return 0.11 if Entrez.api_key else 0.34


def _call(fn: Callable, parse: bool, retries: int = 5, **kw):
    from Bio import Entrez
    socket.setdefaulttimeout(TIMEOUT)   # Bio.Entrez opens URLs without a timeout; a stalled request would hang
    for i in range(retries):
        try:
            h = fn(**kw)
            out = Entrez.read(h) if parse else h.read()
            h.close()
            return out
        except Exception as e:  # noqa: BLE001 - urllib/http/parse errors
            if i == retries - 1:
                raise
            wait = 2 ** i * 2
            print(f"   retry {i + 1}/{retries - 1} ({type(e).__name__}: {str(e)[:80]}) in {wait}s", file=sys.stderr)
            time.sleep(wait)


def resolve_taxon(taxon: str) -> tuple[str, str, str, str]:
    """Name or taxid -> (taxid, scientific name, rank, lineage)."""
    from Bio import Entrez
    if str(taxon).isdigit():
        tid = str(taxon)
    else:
        r = _call(Entrez.esearch, True, db="taxonomy", term=f'"{taxon}"[Scientific Name]')
        if not r["IdList"]:
            r = _call(Entrez.esearch, True, db="taxonomy", term=taxon)
        if not r["IdList"]:
            raise ValueError(f"Taxon '{taxon}' not found in NCBI Taxonomy")
        tid = r["IdList"][0]
    info = _call(Entrez.efetch, True, db="taxonomy", id=tid, retmode="xml")[0]
    return tid, info["ScientificName"], info.get("Rank", ""), info.get("Lineage", "")


def build_query(taxid: str, o: DownloadOptions) -> tuple[str, str]:
    """-> (Entrez query, tag for the output sub-directory)."""
    q = [f"txid{taxid}[Organism:exp]"]
    tag = []
    if not o.all:
        if not o.include_wgs:
            q.append("NOT wgs[filter]")
        if not o.include_mrna:
            q.append("NOT biomol_mrna[PROP]")
        q.append("NOT srcdb_refseq_model[PROP]")
        if not o.include_refseq:
            q.append("NOT refseq[filter]")
        tag += [f"incl-{x}" for x, on in (("wgs", o.include_wgs), ("mrna", o.include_mrna),
                                          ("refseq", o.include_refseq)) if on]
    if o.mitogenome:
        q.append(MITOGENOME)
        tag.append("mitogenome")
    elif o.mito:
        q.append("AND mitochondrion[filter]")
        tag.append("mito")
    if o.gene:
        genes = [canon(g) for g in o.gene.split(",") if g.strip()]
        q.append("AND (" + " OR ".join(clause(g) for g in genes) + ")")
        tag.append("gene-" + "_".join(genes).replace(".", ""))
    if o.minlen or o.maxlen:
        q.append(f"AND {o.minlen or 1}:{o.maxlen or 999999999}[SLEN]")
        tag.append(f"len{o.minlen or 1}-{o.maxlen or 'max'}")
    elif not o.all and not o.mitogenome and not o.include_large:
        q.append(f"AND (1:{LARGE_BP}[SLEN] OR {ORGANELLE})")
    if o.include_large and not o.all:
        tag.append("incl-large")
    if o.query:
        q.append(f"AND ({o.query})")
        tag.append("custom-" + hashlib.sha1(o.query.encode("utf-8")).hexdigest()[:6])
    if o.all:
        tag.insert(0, "all")
    t = o.tag or "-".join(tag) or "markers"
    return " ".join(q), t


def count(term: str) -> int:
    from Bio import Entrez
    return int(_call(Entrez.esearch, True, db="nuccore", term=term, retmax=0)["Count"])


def _n_records(text: str) -> int:
    return text.count("\nLOCUS ") + text.startswith("LOCUS")


_VERSION = re.compile(r"^VERSION\s+(\S+)", re.M)


def _batch_complete(text: str, ids: list[str]) -> bool:
    """A batch file is complete only if it holds exactly the requested accession.versions, each record closed by //."""
    found = _VERSION.findall(text)
    ends = len(re.findall(r"^//\s*$", text, re.M))
    want = {i.strip() for i in ids}
    if all("." in i for i in want):
        return sorted(found) == sorted(want) and ends == len(ids)
    return {f.split(".")[0] for f in found} == {i.split(".")[0] for i in want} and ends == len(ids)


def _efetch_text(ids: list[str]) -> str:
    """One efetch POST with gzip transfer and a timeout (Bio.Entrez does neither)."""
    from Bio import Entrez
    params = {"db": "nuccore", "id": ",".join(ids), "rettype": "gbwithparts", "retmode": "text", "tool": "g2t"}
    if Entrez.email:
        params["email"] = Entrez.email
    if Entrez.api_key:
        params["api_key"] = Entrez.api_key
    req = urllib.request.Request(EFETCH_URL, data=urllib.parse.urlencode(params).encode(),
                                 headers={"Accept-Encoding": "gzip"})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:  # noqa: S310 - fixed NCBI https URL
        raw: bytes = r.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", errors="replace")


def _fetch_batch(ids: list[str]) -> tuple[list[str], str | None]:
    """Fetch one batch; returns (ids, text or None after all attempts failed)."""
    for attempt in range(4):
        try:
            txt = _efetch_text(ids)
            if _batch_complete(txt, ids):
                return ids, txt
            print(f"   batch of {len(ids)}: got {_n_records(txt)} records, retrying", file=sys.stderr)
        except Exception as e:  # noqa: BLE001
            print(f"   batch of {len(ids)}: {type(e).__name__}, retrying", file=sys.stderr)
        time.sleep(2 ** attempt)
    return ids, None


def _previous_list(out: Path) -> list[str]:
    f = out / "accessions.tsv"
    if not f.exists():
        return []
    return [x.strip() for x in f.read_text(encoding="utf-8").splitlines()[1:] if x.strip()]


def _changes(old: list[str], new: list[str]) -> list[tuple[str, str, str]]:
    """(status, accession.version, previous version) for new / updated / removed records."""
    if not old:
        return []
    ob = {a.split(".")[0]: a for a in old}
    nb = {a.split(".")[0]: a for a in new}
    rows = []
    for b, a in nb.items():
        if b not in ob:
            rows.append(("new", a, ""))
        elif ob[b] != a:
            rows.append(("updated", a, ob[b]))
    rows += [("removed", a, "") for b, a in ob.items() if b not in nb]
    return rows


def _acc_list(term: str, total: int, delay: float, workers: int) -> list[str]:
    """accession.version list of a query; pages are fetched in parallel."""
    from Bio import Entrez
    if total == 0:
        return []
    r = _call(Entrez.esearch, True, db="nuccore", term=term, usehistory="y", retmax=0)
    pages: dict[int, list[str]] = {}

    def _page(start: int) -> tuple[int, list[str]]:
        txt = _call(Entrez.efetch, False, db="nuccore", rettype="acc", retmode="text", retstart=start,
                    retmax=5000, webenv=r["WebEnv"], query_key=r["QueryKey"])
        return start, [x.strip() for x in txt.split() if x.strip()]

    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        pfuts = []
        for start in range(0, total, 5000):
            pfuts.append(ex.submit(_page, start))
            time.sleep(delay)
        for pf in as_completed(pfuts):
            pstart, page = pf.result()
            pages[pstart] = page
    return list(dict.fromkeys(a for p in sorted(pages) for a in pages[p]))


SINCE_MARGIN_DAYS = 3   # NCBI may index records a few days after their modification date


def _since_date(since: str, prev: dict, prev_list: list[str], query_changed: bool) -> str:
    """'' (full list) or an Entrez date YYYY/MM/DD for an incremental --since run."""
    if not since:
        return ""
    if not prev_list or query_changed or not prev.get("complete"):
        print("   --since: no complete earlier download of this selection, downloading the full list",
              file=sys.stderr)
        return ""
    if since == "auto":
        d = date.fromisoformat(prev["date"]) - timedelta(days=SINCE_MARGIN_DAYS)
        return d.strftime("%Y/%m/%d")
    return date.fromisoformat(since.replace("/", "-")).strftime("%Y/%m/%d")


def download(taxon: str, out_root: str, options: DownloadOptions | None = None, dry_run: bool = False,
             email: str | None = None, api_key: str | None = None, progress: bool = True,
             store: str | None = None, report: bool = False) -> DownloadResult:
    """Resolve taxon, report composition, and download the selection into <out_root>/<tag>/.

    Records are kept in an accession-level store (``store``, default ``<out_root>/_records.sqlite``);
    only accession.versions missing from it are fetched. ``<out_root>/<tag>/batch_NNNN.gb`` is then
    written from the store, so downstream steps see the usual batch files.
    """
    o = options or DownloadOptions()
    delay = setup_entrez(email, api_key)
    tid, name, rank, lineage = resolve_taxon(taxon)
    term, tag = build_query(tid, o)
    total = count(term)
    base = f"txid{tid}[Organism:exp]"
    comp_q = [("all", ""), ("wgs", " AND wgs[filter]"), ("mrna", " AND biomol_mrna[PROP]"),
              ("refseq_model", " AND srcdb_refseq_model[PROP]"), ("mito", " AND mitochondrion[filter]"),
              ("large", f" AND {LARGE_BP + 1}:999999999[SLEN] NOT {ORGANELLE}")]
    gene_q = [(g, f"({term}) AND {clause(g)}") for g in MITO_GENES + NUCLEAR_GENES]
    # composition and per-gene counts are ~23 slow searches: only for dry runs or report=True
    jobs = ([(k, base + q) for k, q in comp_q] + gene_q) if (dry_run or report) else []
    counts: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=max(1, o.workers)) as ex:
        cfuts = {}
        for k, q in jobs:
            cfuts[ex.submit(count, q)] = k
            time.sleep(delay)
        for cf in as_completed(cfuts):
            counts[cfuts[cf]] = cf.result()
    comp = {k: counts[k] for k, _ in comp_q if k in counts}
    per_gene = {g: counts[g] for g, _ in gene_q if counts.get(g)}
    out = Path(out_root) / tag
    res = DownloadResult(taxid=tid, name=name, rank=rank, tag=tag, query=term, count=total, out_dir=str(out),
                         composition=comp, per_gene=per_gene)
    if total == 0:
        return res
    mfile = out / "manifest.json"
    old_q = json.loads(mfile.read_text(encoding="utf-8")).get("query") if mfile.exists() else None
    res.query_changed = bool(old_q and old_q != term)
    if o.tag and res.query_changed:   # a user-named directory is never silently re-purposed
        raise ValueError(f"{out} holds records of a different query:\n  {old_q}\n"
                         f"current query:\n  {term}\nUse another --tag or output directory.")
    prev = json.loads(mfile.read_text(encoding="utf-8")) if mfile.exists() else {}
    prev_list = _previous_list(out)
    since = _since_date(o.since, prev, prev_list, res.query_changed)
    if since:
        # only records created or modified since the last download (NCBI modification date, [MDAT])
        sterm = f"({term}) AND {since}:3000[MDAT]"
        cand = _acc_list(sterm, count(sterm), delay, o.workers)
        base_new = {a.split(".")[0]: a for a in cand}
        accs = [base_new.pop(a.split(".")[0], a) for a in prev_list] + list(base_new.values())
        res.since = since
        res.candidates = len(cand)
    else:
        # the accession list is always refreshed: NCBI adds and removes records between runs
        accs = _acc_list(term, total, delay, o.workers)
    res.count = len(accs)
    spath = Path(store) if store else Path(out_root) / "_records.sqlite"
    if o.no_store:
        st: RecordStore | DirStore = DirStore(out)
    elif dry_run and not spath.exists():
        res.to_fetch = len(accs)
        return res
    else:
        st = RecordStore(spath)
    missing = st.missing(accs)
    if missing and not dry_run and isinstance(st, RecordStore):
        # records fetched by older g2t versions (plain batch files, any selection) are imported, not re-fetched
        legacy = sorted(Path(out_root).glob("*/batch_*.gb"))
        if legacy:
            st.import_files(legacy, set(missing))
            missing = st.missing(accs)
    res.reused = len(accs) - len(missing)
    res.to_fetch = len(missing)
    if dry_run:
        st.close()
        return res
    out.mkdir(parents=True, exist_ok=True)
    batches = [missing[i:i + o.batch] for i in range(0, len(missing), o.batch)]
    res.batches = len(batches)
    t0 = time.time()
    failed_ids: list[str] = []
    workers = max(1, min(o.workers, len(batches) or 1))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = []
        for ids in batches:
            futs.append(ex.submit(_fetch_batch, ids))
            time.sleep(delay)
        retry: list[list[str]] = []
        for i, fu in enumerate(as_completed(futs), 1):
            ids, txt = fu.result()
            if txt is None:
                retry.append(ids)
            else:
                res.fetched += st.put_many(iter_records(txt.splitlines(keepends=True)))
            if progress:
                print(f"   [{i}/{len(batches)}] {'OK' if txt else 'FAILED, will split'}  {time.time() - t0:.0f}s",
                      file=sys.stderr, flush=True)
    # split failing batches until the failing records are isolated
    while retry:
        ids = retry.pop()
        if len(ids) == 1:
            ids, txt = _fetch_batch(ids)
            if txt is None:
                failed_ids += ids
                res.failed += 1
                print(f"   record {ids[0]}: FAILED", file=sys.stderr)
            else:
                res.fetched += st.put_many(iter_records(txt.splitlines(keepends=True)))
            continue
        for half in (ids[:len(ids) // 2], ids[len(ids) // 2:]):
            half, txt = _fetch_batch(half)
            if txt is None:
                retry.append(half)
            else:
                res.fetched += st.put_many(iter_records(txt.splitlines(keepends=True)))
    failed_set = set(failed_ids)
    have = [a for a in accs if a not in failed_set]
    if isinstance(st, DirStore):
        # store-free: new records were written as new batch files; drop superseded / withdrawn ones
        st.prune(set(have))
    else:
        # write the selection view from the store, in NCBI list order
        nb = (len(have) + o.batch - 1) // o.batch
        unchanged = (not failed_ids and prev.get("complete") and prev.get("query") == term
                     and prev.get("options", {}).get("batch") == o.batch
                     and prev.get("accessions_sha1") == hashlib.sha1("\n".join(accs).encode()).hexdigest()
                     and len(list(out.glob("batch_*.gb"))) == nb)
        for b in range(0 if unchanged else nb):     # nothing new: keep the existing batch files
            f = out / f"batch_{b + 1:04d}.gb"
            tmp = f.with_suffix(".tmp")
            with open(tmp, "w", encoding="utf-8") as fh:
                for a in have[b * o.batch:(b + 1) * o.batch]:
                    fh.write(st.get(a) or "")
            tmp.replace(f)
        for f in sorted(out.glob("batch_*.gb")):
            m = re.fullmatch(r"batch_(\d+)\.gb", f.name)
            if m and int(m.group(1)) > nb:
                f.unlink()
    st.close()
    ch = _changes(prev_list, accs)
    res.new = sum(c[0] == "new" for c in ch)
    res.updated = sum(c[0] == "updated" for c in ch)
    res.removed = sum(c[0] == "removed" for c in ch)
    if ch:
        if res.query_changed:
            ch = [(f"{st_}(query changed)", a, p) for st_, a, p in ch]
        (out / "changes.tsv").write_text("status\taccession\tprevious\n" +
                                         "".join(f"{s}\t{a}\t{p}\n" for s, a, p in ch), encoding="utf-8")
        hist = out / "changes_history.tsv"
        new_h = not hist.exists()
        with open(hist, "a", encoding="utf-8") as fh:
            if new_h:
                fh.write("date\tstatus\taccession\tprevious\n")
            fh.writelines(f"{date.today().isoformat()}\t{s}\t{a}\t{p}\n" for s, a, p in ch)
    elif (out / "changes.tsv").exists():
        (out / "changes.tsv").unlink()
    (out / "accessions.tsv").write_text("accession\n" + "\n".join(accs) + "\n", encoding="utf-8")
    # date of the last check against the full NCBI list (a --since run cannot see withdrawn records)
    full_check = prev.get("last_full_check") if since else date.today().isoformat()
    manifest = dict(taxon=taxon, ncbi_taxid=tid, ncbi_name=name, ncbi_rank=rank, ncbi_lineage=lineage,
                    query=term, tag=tag, date=date.today().isoformat(), count=len(accs),
                    accessions_sha1=hashlib.sha1("\n".join(accs).encode()).hexdigest(),
                    complete=res.failed == 0, failed_batches=res.failed, failed_records=len(failed_ids),
                    fetched=res.fetched, reused=res.reused, new=res.new, updated=res.updated,
                    removed=res.removed, query_changed=res.query_changed,
                    store=None if o.no_store else str(spath), since=since or None,
                    since_candidates=res.candidates if since else None, last_full_check=full_check,
                    taxon_composition=comp, per_gene=per_gene, options=asdict(o))
    mfile.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return res


def describe(res: DownloadResult, dry_run: bool = False) -> list[str]:
    head = f"{res.name} (txid{res.taxid}, {res.rank}) selection '{res.tag}': {res.count} records"
    if dry_run:
        head += f" (dry run: stored {res.reused}, to fetch {res.to_fetch})"
    else:
        head += f" (reused {res.reused}, fetched {res.fetched}, failed {res.failed})"
        head += f"; since last run: new {res.new}, updated {res.updated}, removed {res.removed}"
        if res.query_changed:
            head += " (query changed)"
        if res.since:
            head += f"; modified in NCBI since {res.since}: {res.candidates}"
    lines = [head]
    c = res.composition
    if c:
        lines.append(f"taxon total {c.get('all')}: WGS {c.get('wgs')}, mRNA {c.get('mrna')}, "
                     f"RefSeq models {c.get('refseq_model')}, mitochondrial {c.get('mito')}, "
                     f"nuclear > {LARGE_BP // 1000} kb {c.get('large')}")
    if res.per_gene:
        lines.append("per gene: " + ", ".join(f"{g} {n}" for g, n in res.per_gene.items()))
    return lines + [f"query: {res.query}", f"output: {res.out_dir}"]


def _since_arg(v: str) -> str:
    if v in ("", "auto"):
        return v
    try:
        return date.fromisoformat(v.replace("/", "-")).isoformat()
    except ValueError:
        raise argparse.ArgumentTypeError("expected 'auto' or a date YYYY-MM-DD") from None


def _positive(v: str) -> int:
    n = int(v)
    if n < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return n


def add_selection_args(p: argparse.ArgumentParser) -> None:
    """Selection, update and run options (shared with wrappers that supply taxon/output/credentials)."""
    g = p.add_argument_group("selection (default: marker records, see README)")
    g.add_argument("--all", action="store_true", help="all records, incl. WGS, mRNA and RefSeq")
    g.add_argument("--mito", action="store_true", help="mitochondrial records only")
    g.add_argument("--mitogenome", action="store_true", help="complete mitochondrial genomes only (>= 10 kb)")
    g.add_argument("--gene", default="", metavar="LIST", help="genes, e.g. COI,18S,28S")
    g.add_argument("--minlen", type=_positive, metavar="BP")
    g.add_argument("--maxlen", type=_positive, metavar="BP")
    g.add_argument("--query", default="", metavar="CLAUSE", help="extra Entrez clause, e.g. 'Russia[Country]'")
    g.add_argument("--include-wgs", action="store_true", help="keep WGS contigs")
    g.add_argument("--include-mrna", action="store_true", help="keep mRNA records")
    g.add_argument("--include-refseq", action="store_true", help="keep RefSeq records")
    g.add_argument("--include-large", action="store_true",
                   help=f"keep nuclear records > {LARGE_BP // 1000} kb (organelle records are always kept)")
    g.add_argument("--tag", default="", metavar="NAME", help="output sub-directory (default: from the selection)")

    u = p.add_argument_group("update")
    u.add_argument("--since", type=_since_arg, default="", metavar="auto|DATE",
                   help="fetch only records modified since the last run or DATE (YYYY-MM-DD)")
    x = u.add_mutually_exclusive_group()
    x.add_argument("--store", default=None, metavar="FILE", help="record store (default: OUTPUT/_records.sqlite)")
    x.add_argument("--no-store", action="store_true", help="keep records only in the batch files")

    r = p.add_argument_group("run")
    r.add_argument("-n", "--dry-run", action="store_true", help="report counts and records to fetch; write nothing")
    r.add_argument("--report", action="store_true", help="also count records per gene (extra NCBI searches)")
    r.add_argument("-b", "--batch-size", type=_positive, default=500, metavar="N",
                   help="records per request and batch file (default: %(default)s)")
    r.add_argument("-w", "--workers", type=_positive, default=3, metavar="N",
                   help="parallel NCBI requests (default: %(default)s)")
    r.add_argument("-q", "--quiet", action="store_true", help="no progress output")


def options_from_args(a: argparse.Namespace) -> DownloadOptions:
    return DownloadOptions(all=a.all, mito=a.mito, mitogenome=a.mitogenome, gene=a.gene, minlen=a.minlen,
                           maxlen=a.maxlen, query=a.query, include_wgs=a.include_wgs, include_mrna=a.include_mrna,
                           include_refseq=a.include_refseq, include_large=a.include_large,
                           since=a.since, no_store=a.no_store, batch=a.batch_size, workers=a.workers, tag=a.tag)


EPILOG = """examples:
  g2t-download -t Polynoidae -o gb                 first download (markers)
  g2t-download -t Polynoidae -o gb --since auto    later: only new or modified records
  g2t-download -t Polynoidae -o gb -n              counts only
  g2t-download -t 46593 -o gb --mitogenome

Re-running skips records already downloaded. Exit status 1 if records failed."""


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="g2t-download", description="Download GenBank records of a taxon from NCBI.",
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EPILOG)
    p.add_argument("-t", "--taxon", required=True, metavar="NAME|TAXID", help="taxon name or NCBI taxid")
    p.add_argument("-o", "--output", metavar="DIR", help="output root; records go to DIR/TAG/ (default: <taxon>_gb)")
    p.add_argument("-e", "--email", metavar="ADDR", help="e-mail for NCBI (default: $NCBI_EMAIL)")
    p.add_argument("-k", "--api-key", metavar="KEY", help="NCBI API key (default: $NCBI_API_KEY)")
    p.add_argument("--validate", action="store_true", help="parse the batch files with Biopython afterwards")
    p.add_argument("--resume", action="store_true", help=argparse.SUPPRESS)   # resuming is automatic
    add_selection_args(p)
    a = p.parse_args(argv)
    out = a.output or f"{str(a.taxon).replace(' ', '_')}_gb"
    res = download(a.taxon, out, options_from_args(a), dry_run=a.dry_run, email=a.email, api_key=a.api_key,
                   progress=not a.quiet, store=a.store, report=a.report)
    for line in describe(res, a.dry_run):
        print(line)
    if a.validate and not a.dry_run:
        from Bio import SeqIO
        n = sum(1 for f in sorted(Path(res.out_dir).glob("batch_*.gb")) for _ in SeqIO.parse(str(f), "genbank"))
        print(f"validated: {n} records parse with Biopython")
    return 1 if res.failed else 0


if __name__ == "__main__":
    sys.exit(main())
