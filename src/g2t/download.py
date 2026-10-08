#!/usr/bin/env python3
"""
GenBank downloader (Step 0)
===========================
Download GenBank records for a taxon from NCBI nuccore, ready for ``g2t``.

Compared with the original scripts/download_entrez.py:
  - Default "markers" selection excludes WGS contigs, mRNA, and RefSeq (predicted
    models and NC_/NR_ copies of INSDC records). In many taxa these make up >95%
    of nuccore and are useless for a specimen x gene matrix. ``--all`` keeps them.
  - Selections: ``--mito`` (mitochondrial records), ``--mitogenome`` (complete
    mitochondrial genomes, 10-30 kb), ``--gene COI,18S,28S`` (see g2t.ncbi_genes),
    ``--minlen/--maxlen``, ``--query`` (any extra Entrez clause).
  - Pre-flight report: taxon composition (WGS / mRNA / RefSeq / mito) and number
    of records per gene for the selection; ``--dry-run`` stops there.
  - Downloads by an explicit accession list (WebEnv sessions expire), checks the
    record count of every batch, retries with exponential back-off and resumes
    automatically: re-running the same command skips finished batches.
  - Each selection gets its own sub-directory (``<out>/<tag>/``) with
    accessions.tsv and manifest.json (query, date, counts) for reproducibility.

CLI:
    g2t-download -t Priapulidae -o gb_out                     # markers
    g2t-download -t Priapulidae -o gb_out --mitogenome
    g2t-download -t Priapulidae -o gb_out --gene COI --minlen 500
    g2t-download -t 37891 -o gb_out --dry-run
Email / API key: -e/-k, or environment variables NCBI_EMAIL / NCBI_API_KEY (optional;
with an API key NCBI allows 10 requests/s instead of 3).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import warnings
from dataclasses import asdict, dataclass, field
from datetime import date
from pathlib import Path
from typing import Callable

from g2t.ncbi_genes import MITO_GENES, NUCLEAR_GENES, canon, clause

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
    batch: int = 200
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
    batches: int = 0
    downloaded: int = 0
    skipped: int = 0
    failed: int = 0


def setup_entrez(email: str | None = None, api_key: str | None = None) -> float:
    """Configure Bio.Entrez; returns the delay between requests (s)."""
    from Bio import Entrez
    Entrez.email = email or os.environ.get("NCBI_EMAIL") or None  # type: ignore[assignment]
    Entrez.api_key = api_key or os.environ.get("NCBI_API_KEY") or None  # type: ignore[assignment]
    Entrez.tool = "g2t"
    return 0.11 if Entrez.api_key else 0.34


def _call(fn: Callable, parse: bool, retries: int = 5, **kw):
    from Bio import Entrez
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
        q.append('AND mitochondrion[filter] AND (complete genome[Title] OR "complete mitochondrial genome"[Title] '
                 'OR mitogenome[Title]) AND 10000:30000[SLEN]')
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


def download(taxon: str, out_root: str, options: DownloadOptions | None = None, dry_run: bool = False,
             email: str | None = None, api_key: str | None = None, progress: bool = True) -> DownloadResult:
    """Resolve taxon, report composition, and download the selection into <out_root>/<tag>/."""
    from Bio import Entrez
    o = options or DownloadOptions()
    delay = setup_entrez(email, api_key)
    tid, name, rank, lineage = resolve_taxon(taxon)
    term, tag = build_query(tid, o)
    total = count(term)
    base = f"txid{tid}[Organism:exp]"
    comp = {k: count(base + q) for k, q in [("all", ""), ("wgs", " AND wgs[filter]"),
                                             ("mrna", " AND biomol_mrna[PROP]"),
                                             ("refseq_model", " AND srcdb_refseq_model[PROP]"),
                                             ("mito", " AND mitochondrion[filter]")]}
    per_gene = {}
    for g in MITO_GENES + NUCLEAR_GENES:
        n = count(f"({term}) AND {clause(g)}")
        if n:
            per_gene[g] = n
        time.sleep(delay)
    out = Path(out_root) / tag
    res = DownloadResult(taxid=tid, name=name, rank=rank, tag=tag, query=term, count=total, out_dir=str(out),
                         composition=comp, per_gene=per_gene)
    if dry_run or total == 0:          # nothing is written for a dry run
        return res
    mfile = out / "manifest.json"
    if mfile.exists():
        old_q = json.loads(mfile.read_text(encoding="utf-8")).get("query")
        if old_q and old_q != term:
            raise ValueError(f"{out} holds records of a different query:\n  {old_q}\n"
                             f"current query:\n  {term}\nUse another --tag or output directory.")
    out.mkdir(parents=True, exist_ok=True)
    # the accession list is always refreshed: NCBI adds and removes records between runs
    r = _call(Entrez.esearch, True, db="nuccore", term=term, usehistory="y", retmax=0)
    accs = []
    for start in range(0, total, 5000):
        txt = _call(Entrez.efetch, False, db="nuccore", rettype="acc", retmode="text", retstart=start,
                    retmax=5000, webenv=r["WebEnv"], query_key=r["QueryKey"])
        accs += [x.strip() for x in txt.split() if x.strip()]
        time.sleep(delay)
    accs = list(dict.fromkeys(accs))
    res.count = len(accs)
    (out / "accessions.tsv").write_text("accession\n" + "\n".join(accs) + "\n", encoding="utf-8")
    nb = (len(accs) + o.batch - 1) // o.batch
    res.batches = nb
    t0 = time.time()
    for b in range(nb):
        ids = accs[b * o.batch:(b + 1) * o.batch]
        f = out / f"batch_{b + 1:04d}.gb"
        if f.exists() and _batch_complete(f.read_text(errors="ignore"), ids):
            res.skipped += 1
            continue
        ok = False
        for attempt in range(4):
            try:
                txt = _call(Entrez.efetch, False, db="nuccore", id=",".join(ids), rettype="gbwithparts", retmode="text")
                if _batch_complete(txt, ids):
                    f.write_text(txt, encoding="utf-8")
                    ok = True
                    break
                print(f"   batch {b + 1}: got {_n_records(txt)}/{len(ids)} records, retrying", file=sys.stderr)
            except Exception as e:  # noqa: BLE001
                print(f"   batch {b + 1}: {type(e).__name__}, retrying", file=sys.stderr)
            time.sleep(2 ** attempt * 2)
        res.downloaded += ok
        res.failed += not ok
        if progress:
            print(f"   [{b + 1}/{nb}] {'OK' if ok else 'FAILED'}  {time.time() - t0:.0f}s", flush=True)
        time.sleep(delay)
    # batch files beyond the current list (smaller list or larger batch size) would duplicate records downstream
    for f in sorted(out.glob("batch_*.gb")):
        m = re.fullmatch(r"batch_(\d+)\.gb", f.name)
        if m and int(m.group(1)) > nb:
            f.unlink()
    manifest = dict(taxon=taxon, ncbi_taxid=tid, ncbi_name=name, ncbi_rank=rank, ncbi_lineage=lineage,
                    query=term, tag=tag, date=date.today().isoformat(), count=len(accs),
                    accessions_sha1=hashlib.sha1("\n".join(accs).encode()).hexdigest(),
                    complete=res.failed == 0, failed_batches=res.failed,
                    taxon_composition=comp, per_gene=per_gene, options=asdict(o))
    mfile.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return res


def describe(res: DownloadResult, dry_run: bool = False) -> list[str]:
    c = res.composition
    genes = ", ".join(f"{g} {n}" for g, n in res.per_gene.items()) or "none"
    head = f"{res.name} (txid{res.taxid}, {res.rank}) selection '{res.tag}': {res.count} records"
    if dry_run:
        head += " (dry run, nothing downloaded)"
    else:
        head += f", {res.batches} batches (new {res.downloaded}, existing {res.skipped}, failed {res.failed})"
    return [head,
            f"taxon total {c.get('all')}: WGS {c.get('wgs')}, mRNA {c.get('mrna')}, "
            f"RefSeq models {c.get('refseq_model')}, mitochondrial {c.get('mito')}",
            f"per gene: {genes}",
            f"query: {res.query}",
            f"output: {res.out_dir}"]


def add_selection_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("selection")
    g.add_argument("--all", action="store_true", help="All records incl. WGS/mRNA/RefSeq (can be huge)")
    g.add_argument("--mito", action="store_true", help="Mitochondrial records only")
    g.add_argument("--mitogenome", action="store_true", help="Complete mitochondrial genomes only (10-30 kb)")
    g.add_argument("--gene", default="", help="Comma-separated genes, e.g. COI,18S,28S (see g2t/ncbi_genes.py)")
    g.add_argument("--minlen", type=int)
    g.add_argument("--maxlen", type=int)
    g.add_argument("--query", default="", help="Extra Entrez clause, e.g. 'Russia[Country]'")
    g.add_argument("--include-wgs", action="store_true")
    g.add_argument("--include-mrna", action="store_true")
    g.add_argument("--include-refseq", action="store_true")
    g.add_argument("--tag", default="", help="Name of the output sub-directory (default: from the selection)")
    g.add_argument("--dry-run", action="store_true", help="Only report counts and query")
    g.add_argument("-b", "--batch-size", type=int, default=200, help="Records per batch file (default 200)")


def options_from_args(a: argparse.Namespace) -> DownloadOptions:
    return DownloadOptions(all=a.all, mito=a.mito, mitogenome=a.mitogenome, gene=a.gene, minlen=a.minlen,
                           maxlen=a.maxlen, query=a.query, include_wgs=a.include_wgs, include_mrna=a.include_mrna,
                           include_refseq=a.include_refseq, batch=a.batch_size, tag=a.tag)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Download GenBank records of a taxon from NCBI (markers by default)",
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("-t", "--taxon", required=True, help="Taxon name or NCBI taxid")
    p.add_argument("-o", "--output", default=None,
                   help="Output root (default: <taxon>_gb/); files go to <output>/<tag>/")
    p.add_argument("-e", "--email", default=None, help="Email for NCBI (or env NCBI_EMAIL)")
    p.add_argument("-k", "--api-key", default=None, help="NCBI API key (or env NCBI_API_KEY)")
    p.add_argument("--resume", action="store_true", help="(kept for compatibility; resuming is automatic)")
    p.add_argument("--validate", action="store_true", help="Parse downloaded files with Biopython afterwards")
    add_selection_args(p)
    a = p.parse_args(argv)
    out = a.output or f"{str(a.taxon).replace(' ', '_')}_gb"
    res = download(a.taxon, out, options_from_args(a), dry_run=a.dry_run, email=a.email, api_key=a.api_key)
    for line in describe(res, a.dry_run):
        print(line)
    if a.validate and not a.dry_run:
        from Bio import SeqIO
        n = sum(1 for f in sorted(Path(res.out_dir).glob("batch_*.gb")) for _ in SeqIO.parse(str(f), "genbank"))
        print(f"validated: {n} records parse with Biopython")
    return 1 if res.failed else 0


if __name__ == "__main__":
    sys.exit(main())
