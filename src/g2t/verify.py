"""g2t-verify: check gene assignments against the sequence itself (rRNA profile HMMs).

g2t assigns gene types from the record title (DEFINITION). Titles are sometimes wrong — e.g. five
Hesionidae 16S records (MG523356-60) are titled "small subunit ribosomal RNA ...; mitochondrial", so
they land in the 12s column — and multi-region rDNA records ("18S partial; ITS1, 5.8S, ITS2; 28S
partial") are put in the 18s and 28s columns although they hold only primer-end fragments of 18S/28S
(Enchytraeidae: at most 67 bp). See docs/evaluation/rdna-regions.md.

g2t-verify scans the sequences with nucleotide profile HMMs of mitochondrial 12S and 16S and nuclear
18S, 5.8S and 28S rRNA (HMMER nhmmer via pyhmmer), locates each region, derives ITS1/ITS2 from the
5.8S and its flanks, and reports per record which gene types the sequence supports:

    agree       the regions found give the same gene types as g2t
    partly      some of g2t's types are supported, others not (e.g. 28s,its1-its2 -> its1-its2)
    different   the sequence is another rRNA than g2t says (e.g. 12s -> 16s)
    no_hit      g2t says rRNA, but no rRNA model matches (check the record)
    unverifiable  spacer only (e.g. "ITS2 and 28S, partial", 237 bp): ITS has no model, and without
                a 5.8S there is nothing to place it by (a short 18S / 28S end is not enough);
                g2t's type is kept
    new         g2t found no gene type, the sequence is an rRNA

A region counts towards a gene type when it is >= --min-bp long (default 150), covers >= --min-cover
of the model (default 0.5), or makes up at least half of a short record (>= 50 bp); ITS1 / ITS2 count
from 50 bp. Overlapping hits of different models (mitochondrial 16S also matches nuclear LSU)
compete: the higher score wins.

Speed (Enchytraeidae, 4,750 records: 11 min 36 s -> 1 min 23 s (CPU 938 s -> 162 s) on a busy 16-core
machine, 9 s when rerun; verdicts the same or better):
    - identical sequences are scanned once; (model, chunk) jobs run on all cores
    - models longer than 600 nodes (12S, 16S, 18S, 28S) are scanned as overlapping 300-node
      tiles whose hits are joined; a region end the tiles do not reach is aligned with the full
      model on that stretch only
    - cascade: the models of g2t's region run first, the others only for records they do not
      explain (strong hits >= 0.3 bits/bp covering all but 50 bp); --all-models runs all
    - hits are cached per sequence and model (~/.cache/g2t/verify_hits.sqlite; --no-cache)
    - records without any hit are scanned again at E <= 1e-3; E-values are per 1 Mb of target,
      so results do not depend on the batch

Default models: nuclear 18S, 5.8S, 28S = the profile HMMs inside Rfam 15.1 RF01960, RF00002,
RF02543 (CC0, bundled: data/rrna_rfam15.1.hmm.gz); mitochondrial 12S, 16S = barrnap 0.9 mito.hmm
(GPL-3.0: downloaded on first use, checked by SHA-256, cached; not bundled). --models FILE uses your
own HMM file(s) (HMM names containing 12S, 16S, 18S or SSU, 5_8S or 5.8S, 28S or LSU).

CLI:
    g2t-verify -i gb/markers --g2t out -o verify           # rRNA, unmatched and conflicted records
    g2t-verify -i gb/markers --g2t out -o verify --all     # every record
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import logging
import os
import re
import sys
import urllib.request
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger("g2t.verify")

BARRNAP = "https://raw.githubusercontent.com/tseemann/barrnap/f843d17d39561eda7ce27cc43e3c6a8104aeb9ce/db/"  # tag 0.9
DEFAULT_MODELS = {  # downloaded (GPL-3.0, not bundled)
    "mito.hmm": "cf70583c8d0a54c7343b8426503a107509cb37bb883bdea8901efa1e25ce8e01",
}
BUNDLED = "data/rrna_rfam15.1.hmm.gz"  # Rfam 15.1 RF01960 / RF00002 / RF02543 profile HMMs (CC0)
RRNA_TYPES = {"12s", "16s", "18s", "28s", "its1-its2", "18-28s"}
REGION_ORDER = ["18s", "its1", "5.8s", "its2", "28s", "12s", "16s"]
EVALUE = 1e-5
TILE, TILE_STEP = 300, 200  # long models are scanned as overlapping tiles of 300 nodes
TILE_ABOVE = 600  # nodes: 12S, 16S, 18S, 28S are tiled; 5.8S (154) is not
RESCUE_EVALUE = 1e-3  # records without any hit (Enchytraeidae: 10 ITS records, all weak 5.8S, types as g2t)
Z_MB = 1.0  # E-values per 1 Mb of target, whatever the batch: chunked and cached scans give the same hits
CACHE_VERSION = 1
# cascade: the models of the region(s) g2t assigned run first; the others only when the sequence is
# not explained by them (calibrated on Hesionidae and Enchytraeidae: cross-model hits, e.g. LSU on
# mitochondrial 16S, score <= 0.29 bits/bp; same final types for all 5,654 records, half the work)
STRONG_BITS_PER_BP = 0.3
MAX_UNEXPLAINED_BP = 50
STAGE1 = {"12s": {"12s"}, "16s": {"16s"}, "18s": {"18s"}, "28s": {"28s"},
          "its1-its2": {"18s", "5.8s", "28s"}, "18-28s": {"18s", "5.8s", "28s"}}


def model_region(name: str) -> str | None:
    """HMM name -> region: 12S_rRNA -> 12s, 5_8S_rRNA -> 5.8s; 5S and unknown names -> None."""
    n = name.lower()
    for pat, reg in ((r"5[._]8s", "5.8s"), (r"(^|[^0-9])12s", "12s"), (r"(^|[^0-9])16s", "16s"),
                     (r"(^|[^0-9])18s|ssu", "18s"), (r"(^|[^0-9])2[5-8]s|lsu", "28s")):
        if re.search(pat, n):
            return reg
    return None


# ------------------------------------------------------------------ models


def cache_dir() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "g2t" / "models"


def nuclear_models(directory: Path | None = None) -> Path:
    """The bundled Rfam 15.1 18S / 5.8S / 28S HMMs, unpacked once into ``directory`` (default: cache)."""
    import gzip
    from importlib import resources

    d = directory or cache_dir()
    d.mkdir(parents=True, exist_ok=True)
    nuc = d / "rrna_rfam15.1.hmm"
    if not nuc.exists():
        nuc.write_bytes(gzip.decompress(resources.files("g2t").joinpath(BUNDLED).read_bytes()))
    return nuc


def default_models(directory: Path | None = None) -> list[Path]:
    """Rfam 15.1 nuclear rRNA HMMs (bundled) + barrnap 0.9 mito.hmm (downloaded once, SHA-256 checked)."""
    d = directory or cache_dir()
    out = [nuclear_models(d)]
    for name, sha in DEFAULT_MODELS.items():
        f = d / f"barrnap-0.9-{name}"
        if not f.exists() or hashlib.sha256(f.read_bytes()).hexdigest() != sha:
            logger.info("downloading %s (barrnap 0.9, GPL-3.0) to %s", name, f)
            req = urllib.request.Request(BARRNAP + name, headers={"User-Agent": "g2t-verify"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            if hashlib.sha256(data).hexdigest() != sha:
                raise RuntimeError(f"{name}: checksum mismatch, not used")
            f.write_bytes(data)
        out.append(f)
    return out


# ------------------------------------------------------------------ scanning


@dataclass
class Hit:
    region: str
    start: int  # 1-based, on the record
    end: int
    hmm_from: int
    hmm_to: int
    hmm_len: int
    score: float
    strand: str = "+"

    @property
    def cover(self) -> float:
        return (self.hmm_to - self.hmm_from + 1) / self.hmm_len if self.hmm_len else 0.0


def _pyhmmer() -> Any:
    try:
        import pyhmmer
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise SystemExit("g2t-verify needs pyhmmer: pip install pyhmmer  (or pip install -e '.[verify]')") from e
    return pyhmmer


@dataclass
class Model:
    region: str
    hmm: Any  # pyhmmer.plan7.HMM
    key: str  # file checksum + HMM name: identifies the model in the hit cache
    offset: int = 0  # a tile: its node 1 is node offset + 1 of the full model
    full: int = 0  # a tile: length of the full model


def load_models(model_files: list[Path]) -> list[Model]:
    from pyhmmer.plan7 import HMMFile

    _pyhmmer()
    out = []
    for f in model_files:
        digest = hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]
        with HMMFile(str(f)) as hf:
            for h in hf:
                region = model_region(_s(h.name))
                if region:
                    out.append(Model(region, h, f"{digest}:{_s(h.name)}"))
    if not out:
        raise SystemExit("no 12S/16S/18S/5.8S/28S models in " + ", ".join(map(str, model_files)))
    return out


class HitCache:
    """Hits per (sequence, model, mode) in SQLite: a rerun scans only sequences it has not seen.

    Hits do not depend on the batch because the E-value database size is fixed (Z_MB)."""

    def __init__(self, path: Path | None):
        self.db = None
        if path is not None:
            import sqlite3

            path.parent.mkdir(parents=True, exist_ok=True)
            self.db = sqlite3.connect(str(path), check_same_thread=False)
            # a cache: losing the last writes on a crash only means rescanning them
            self.db.execute("PRAGMA journal_mode=WAL")
            self.db.execute("PRAGMA synchronous=OFF")
            self.db.execute("CREATE TABLE IF NOT EXISTS hits (k TEXT PRIMARY KEY, v TEXT)")

    @staticmethod
    def key(seq: str, model: Model, evalue: float) -> str:
        return hashlib.sha1(f"{CACHE_VERSION}|{model.key}|{evalue:g}|{seq}".encode()).hexdigest()

    def get(self, keys: list[str]) -> dict[str, list[Hit]]:
        if self.db is None or not keys:
            return {}
        import json

        out = {}
        for i in range(0, len(keys), 500):
            part = keys[i:i + 500]
            q = f"SELECT k, v FROM hits WHERE k IN ({','.join('?' * len(part))})"
            for k, v in self.db.execute(q, part):
                out[k] = [Hit(*x) for x in json.loads(v)]
        return out

    def put(self, items: dict[str, list[Hit]]) -> None:
        if self.db is None or not items:
            return
        import json

        self.db.executemany("INSERT OR REPLACE INTO hits VALUES (?, ?)", [
            (k, json.dumps([[h.region, h.start, h.end, h.hmm_from, h.hmm_to, h.hmm_len, h.score, h.strand]
                            for h in v])) for k, v in items.items()])
        self.db.commit()


def hit_cache_path() -> Path:
    return cache_dir().parent / "verify_hits.sqlite"


def _nhmmer(model: Model, seqs: list[str], evalue: float) -> list[list[Hit]]:
    """One model against a chunk of sequences, single thread; hits per sequence (same order)."""
    pyhmmer = _pyhmmer()
    from pyhmmer.easel import DigitalSequenceBlock, TextSequence

    al = model.hmm.alphabet
    block = DigitalSequenceBlock(al, [TextSequence(name=str(i).encode(), sequence=s).digitize(al)
                                      for i, s in enumerate(seqs)])
    out: list[list[Hit]] = [[] for _ in seqs]
    for hits in pyhmmer.hmmer.nhmmer([model.hmm], block, cpus=1, E=evalue, Z=Z_MB):
        for h in hits:
            if not h.included:
                continue
            for d in h.domains:
                if not d.included:
                    continue
                a, b = sorted((d.env_from, d.env_to))
                out[int(_s(h.name))].append(Hit(
                    model.region, a, b, d.alignment.hmm_from + model.offset, d.alignment.hmm_to + model.offset,
                    model.full or model.hmm.M, d.score, "+" if d.env_from <= d.env_to else "-"))
    return out


def scan(seqs: dict[str, str], model_files: list[Path] | list[Model], cpus: int = 0, evalue: float = EVALUE,
         want: dict[str, set[str]] | None = None, cache: HitCache | None = None) -> dict[str, list[Hit]]:
    """nhmmer of the models against the sequences (pyhmmer); hits with E <= ``evalue`` per Z_MB of target.

    Identical sequences are scanned once; (model, chunk) jobs run on ``cpus`` threads (default: all
    cores; HMMER releases the GIL). ``want``: per record, the regions whose models to run (default:
    all)."""
    from concurrent.futures import ThreadPoolExecutor

    models = [m for m in model_files if isinstance(m, Model)] or load_models(model_files)  # type: ignore[arg-type]
    cache = cache or HitCache(None)
    uniq: dict[str, list[str]] = {}  # cleaned sequence -> records
    for n, s in seqs.items():
        uniq.setdefault(_clean(s).upper(), []).append(n)
    per_seq: dict[str, list[Hit]] = {s: [] for s in uniq}
    todo: list[tuple[Model, list[str]]] = []
    for m in models:
        targets = [s for s, names in uniq.items()
                   if want is None or any(m.region in want.get(n, ()) for n in names)]
        keys = {s: cache.key(s, m, evalue) for s in targets}
        got = cache.get(list(keys.values()))
        for s in targets:
            if keys[s] in got:
                per_seq[s] += got[keys[s]]
        rest = [s for s in targets if keys[s] not in got]
        if rest:
            todo.append((m, rest))
    threads = cpus or os.cpu_count() or 1
    jobs: list[tuple[Model, list[str]]] = []
    work = {id(r): m.hmm.M * sum(map(len, r)) for m, r in todo}  # cost ~ model length x residues
    total = sum(work.values()) or 1
    for m, rest in todo:  # ~4 jobs per thread, split by cost; longest sequences dealt first
        k = max(1, min(len(rest), round(4 * threads * work[id(rest)] / total)))
        rest = sorted(rest, key=len, reverse=True)
        jobs += [(m, rest[i::k]) for i in range(k)]
    jobs.sort(key=lambda j: -j[0].hmm.M * sum(map(len, j[1])))  # big jobs first: no long tail
    if jobs:
        logger.info("nhmmer: %d models x %d sequences, %d jobs on %d threads%s",
                    len(todo), len({s for _, r in todo for s in r}), len(jobs), threads,
                    f" (E <= {evalue:g})" if evalue != EVALUE else "")
    with ThreadPoolExecutor(threads) as ex:
        for (m, chunk), res in zip(jobs, ex.map(lambda j: _nhmmer(j[0], j[1], evalue), jobs)):
            cache.put({cache.key(s, m, evalue): h for s, h in zip(chunk, res)})
            for s, h in zip(chunk, res):
                per_seq[s] += h
    tiled = {m.region for m in models if m.full}
    return {n: merge_tiles(per_seq[s], tiled) for s, names in uniq.items() for n in names if per_seq[s]}


def _s(x: Any) -> str:
    return x.decode() if isinstance(x, bytes) else str(x)


def _clean(s: str) -> str:
    return re.sub(r"[^ACGTUNRYKMSWBDHVacgtunrykmswbdhv]", "N", s)


# ------------------------------------------------------------------ decisions


@dataclass
class Regions:
    spans: dict[str, tuple[int, int, float]] = field(default_factory=dict)  # region -> (start, end, cover)
    strand: str = "+"
    length: int = 0  # of the record


def regions(hits: list[Hit], length: int) -> Regions:
    """Best hit per region (by score); ITS1 / ITS2 from the 5.8S and its 18S / 28S flanks.

    Hits of different models that overlap by more than half of the shorter one compete: the higher
    score wins (mitochondrial 16S also aligns to the nuclear LSU model, 12S to SSU)."""
    r = Regions(length=length)
    kept: list[Hit] = []
    for h in sorted(hits, key=lambda x: -x.score):
        if all(h.region == k.region or _overlap(h, k) <= 0.5 for k in kept):
            kept.append(h)
    best: dict[str, Hit] = {}
    for h in kept:
        if h.region and (h.region not in best or h.score > best[h.region].score):
            best[h.region] = h
    if best:
        r.strand = max(best.values(), key=lambda h: h.score).strand
    for k, h in best.items():
        r.spans[k] = (h.start, h.end, round(h.cover, 3))
    if "5.8s" in best:
        s58, e58 = best["5.8s"].start, best["5.8s"].end
        a = best["18s"].end + 1 if "18s" in best and best["18s"].end < s58 else 1
        b = best["28s"].start - 1 if "28s" in best and best["28s"].start > e58 else length
        if s58 - a >= 1:
            r.spans["its1"] = (a, s58 - 1, 0.0)
        if b - e58 >= 1:
            r.spans["its2"] = (e58 + 1, b, 0.0)
    return r


def stage1_regions(g2t: list[str]) -> set[str] | None:
    """Models to run first for g2t's types (always 5.8S, it is cheap); None = all models."""
    if not g2t or any(t not in STAGE1 for t in g2t):
        return None
    return {"5.8s"}.union(*(STAGE1[t] for t in g2t))


def explained(hits: list[Hit], length: int) -> bool:
    """Strong winning hits (>= STRONG_BITS_PER_BP), plus ITS placed by a strong 5.8S, cover all but
    MAX_UNEXPLAINED_BP of the record: no other model needs to run."""
    r = regions(hits, length)
    covered = [False] * length
    strong = set()
    for h in hits:
        if h.region in r.spans and r.spans[h.region][:2] == (h.start, h.end) \
                and h.score >= STRONG_BITS_PER_BP * (h.end - h.start + 1):
            strong.add(h.region)
            covered[h.start - 1: h.end] = [True] * (h.end - h.start + 1)
    if "5.8s" in strong:
        for k in ("its1", "its2"):
            if k in r.spans:
                a, b, _ = r.spans[k]
                covered[a - 1: b] = [True] * (b - a + 1)
    return bool(strong) and length - sum(covered) <= MAX_UNEXPLAINED_BP


def ragged(hits: list[Hit], length: int, tiled: set[str], slack: int = 20, overlap: int = 200
           ) -> list[tuple[str, int, int]]:
    """Ends of tiled regions that stop >= ``slack`` bp short of the next region (or the record end)
    on a side where the model has >= ``slack`` nodes left: weakly conserved ends that single tiles
    miss and the full model aligns (KU728933 16S: tiles 563-1206, full model 55-1206).
    Returns (region, from, to): the stretch to align with the full model (gap + ``overlap`` bp)."""
    r = regions(hits, length)
    out = []
    for k, (a, b, _) in r.spans.items():
        if k not in tiled:
            continue
        h = next(x for x in hits if x.region == k and (x.start, x.end) == (a, b))
        left = max([y for x, (_, y, _) in r.spans.items() if x != k and y < a], default=0)
        right = min([x0 for x, (x0, _, _) in r.spans.items() if x != k and x0 > b], default=length + 1)
        five, three = (h.hmm_from - 1, h.hmm_len - h.hmm_to)  # model nodes left on each side
        if h.strand == "-":
            five, three = three, five
        if a - left - 1 >= slack and five >= slack:
            out.append((k, left + 1, min(b, a + overlap)))
        if right - b - 1 >= slack and three >= slack:
            out.append((k, max(a, b - overlap), right - 1))
    return out


def extend(hits: list[Hit], region: str, ext: list[Hit]) -> list[Hit]:
    """Widen the region's best hit by full-model hits (same strand) that overlap it."""
    own = [h for h in hits if h.region == region]
    if not own or not ext:
        return hits
    best = max(own, key=lambda h: h.score)
    new = best
    for e in ext:
        if e.strand == best.strand and min(e.end, new.end) - max(e.start, new.start) >= 20:
            new = Hit(region, min(new.start, e.start), max(new.end, e.end), min(new.hmm_from, e.hmm_from),
                      max(new.hmm_to, e.hmm_to), new.hmm_len, new.score, new.strand)
    return [new if h is best else h for h in hits]


def _overlap(a: Hit, b: Hit) -> float:
    o = min(a.end, b.end) - max(a.start, b.start) + 1
    return max(0, o) / min(a.end - a.start + 1, b.end - b.start + 1)


def supported_types(r: Regions, min_bp: int = 150, min_cover: float = 0.5, min_its: int = 50) -> list[str]:
    """g2t gene types the sequence supports."""
    s = r.spans

    def ok(k: str) -> bool:
        if k not in s:
            return False
        n = s[k][1] - s[k][0] + 1
        # long enough, or most of the model, or most of a short record (a 175 bp 28S fragment is 28S)
        return n >= min_bp or s[k][2] >= min_cover or (n >= 50 and r.length > 0 and n >= 0.5 * r.length)

    def its_ok() -> bool:
        return any(k in s and s[k][1] - s[k][0] + 1 >= min_its for k in ("its1", "its2"))

    types = [k for k in ("12s", "16s", "18s", "28s") if ok(k)]
    if its_ok():
        types.append("its1-its2")
    if "18s" in types and "28s" in types and "its1-its2" in types:
        types = [t for t in types if t not in ("18s", "28s", "its1-its2")] + ["18-28s"]
    return types


def g2t_types(gene_type: Any) -> list[str]:
    if not isinstance(gene_type, str) or not gene_type.strip():
        return []
    return [t.strip() for t in gene_type.split(",") if t.strip()]


def verdict(g2t: list[str], seq: list[str], any_hit: bool) -> str:
    g, s = set(g2t), set(seq)
    for x in (g, s):  # 18-28s = 18s + its1-its2 + 28s
        if "18-28s" in x:
            x |= {"18s", "28s", "its1-its2"}
            x.discard("18-28s")
    rg = g & (RRNA_TYPES | {"18s", "28s", "its1-its2"})
    if not rg:
        return "new" if s else ("agree" if not any_hit else "check")
    if not s:
        if "its1-its2" in rg:  # ITS without a 5.8S anchor (at most a short 18S / 28S end): cannot be checked
            return "unverifiable"
        return "no_hit"
    if rg == s:
        return "agree"
    if rg & s:
        return "partly"
    return "different"


# ------------------------------------------------------------------ input / output


def read_status(g2t_dir: str | Path) -> pd.DataFrame:
    """record_status.csv of a g2t run (every input record with its gene type)."""
    p = Path(g2t_dir)
    for f in (p / "labeled_genes" / "record_status.csv", p / "record_status.csv"):
        if f.exists():
            t = pd.read_csv(f, low_memory=False, dtype=str)
            acc = next(c for c in ("ACCESSION", "accession", "LocusID") if c in t.columns)
            t["acc_base"] = t[acc].str.split(".").str[0]
            return t
    raise SystemExit(f"{g2t_dir}: no labeled_genes/record_status.csv (run g2t first)")


_VERSION = re.compile(r"^VERSION {5}(\S+)", re.M)
_ACCESSION = re.compile(r"^ACCESSION {3}(\S+)", re.M)
_ORGANISM = re.compile(r"^  ORGANISM  (.+)$", re.M)
_NOT_BASES = str.maketrans("", "", "0123456789 \t\r\n/")


def read_sequences(paths: Iterable[str | Path], keep: set[str] | None = None) -> Iterator[tuple[str, str, str]]:
    """(accession.version, organism, sequence) of GenBank records; records without sequence (CONTIG,
    WGS / TSA masters) are skipped. Reads each file whole and cuts it at "//": ~20 times faster than
    Biopython's parser (same output on Enchytraeidae and Hesionidae)."""
    from g2t._extract.files import open_genbank
    from g2t.mitostats import _gb_files

    for f in _gb_files(paths):
        with open_genbank(str(f)) as fh:
            text = fh.read()
        for rec in text.split("\n//"):
            i = rec.find("\nORIGIN")
            if i < 0:
                continue
            m = _VERSION.search(rec, 0, i) or _ACCESSION.search(rec, 0, i)
            if not m or (keep is not None and m.group(1).split(".")[0] not in keep):
                continue
            j = rec.find("\n", i + 1)
            s = rec[j + 1:].translate(_NOT_BASES).upper() if j > 0 else ""
            if s:
                o = _ORGANISM.search(rec, 0, i)
                yield m.group(1), o.group(1).strip() if o else "", s



def select(status: pd.DataFrame, all_records: bool = False) -> set[str]:
    if all_records:
        return set(status["acc_base"])
    gt = status["gene_type"].fillna("")
    rrna = gt.apply(lambda g: bool(set(g2t_types(g)) & (RRNA_TYPES | {"18s", "28s"})))
    other = status["status"].fillna("") != "assigned"
    if "Conflict" in status.columns:
        other |= status["Conflict"].fillna("").str.lower() == "true"
    return set(status.loc[rrna | other, "acc_base"])


def _fmt(r: Regions) -> str:
    order = sorted(r.spans.items(), key=lambda kv: kv[1][0])
    return "; ".join(f"{k}:{a}-{b}" + (f" ({c:.0%} of model)" if c else "") for k, (a, b, c) in order)


def _sub_hmms(model: Model, spans: list[tuple[int, int]]) -> list[Any]:
    """Profile HMMs of node ranges of ``model``, built with pyhmmer's Builder (so E-values are
    calibrated) from 2,000 sequences sampled from the full model with a fixed seed; cached on disk."""
    from pyhmmer.easel import TextMSA, TextSequence
    from pyhmmer.plan7 import Background, Builder, HMMFile

    tag = hashlib.sha1(f"{model.key}:{spans}:v1".encode()).hexdigest()[:16]
    f = cache_dir() / f"tiles-{tag}.hmm"
    if f.exists():
        with HMMFile(str(f)) as hf:
            got = list(hf)
        if len(got) == len(spans):
            return got
    logger.info("building %d tiles of %s (once)", len(spans), _s(model.hmm.name))
    msa = model.hmm.emit_alignment(2000, randomness=42).textize()
    cols = [i for i, c in enumerate(msa.reference) if c not in ".-~"]  # alignment column of each node
    al = model.hmm.alphabet
    builder = Builder(al, architecture="hand", effective_number=2000.0, seed=42)
    out = []
    for lo, hi in spans:
        a, b = cols[lo - 1], cols[hi - 1] + 1
        sub = TextMSA(name=f"{_s(model.hmm.name)}_{lo}-{hi}".encode(), sequences=[
            TextSequence(name=f"s{i}".encode(), sequence=r[a:b]) for i, r in enumerate(msa.alignment)])
        sub.reference = msa.reference[a:b]
        out.append(builder.build_msa(sub.digitize(al), Background(al))[0])
    f.parent.mkdir(parents=True, exist_ok=True)
    with open(f, "wb") as fh:
        for h in out:
            h.write(fh)
    return out


def tile_models(models: list[Model]) -> list[Model]:
    """Models longer than TILE_ABOVE nodes -> overlapping tiles of TILE nodes (step TILE_STEP).

    nhmmer's cost on a real hit grows with model length x hit length; a 3,400 bp 28S record against
    the 3,401-node 28S model is ~0.3 CPU-s. Each tile aligns only around its own part, and its hits
    are joined back into one region (merge_tiles)."""
    out = []
    for m in models:
        M = m.hmm.M
        if M <= TILE_ABOVE:
            out.append(m)
            continue
        starts = list(range(1, M - TILE + 1, TILE_STEP))
        if starts[-1] + TILE - 1 < M:
            starts.append(M - TILE + 1)
        spans = [(lo, lo + TILE - 1) for lo in starts]
        for (lo, hi), h in zip(spans, _sub_hmms(m, spans)):
            out.append(Model(m.region, h, f"{m.key}:{lo}-{hi}:v1", offset=lo - 1, full=M))
    return out


def merge_tiles(hits: list[Hit], tiled: set[str] | frozenset[str] = frozenset({"12s", "16s", "18s", "28s"})
                ) -> list[Hit]:
    """Join hits of one region that continue each other on the record and on the model (tiles).

    Hits are joined across gaps (variable expansion segments of 28S, a weak tile) when the record
    gap is at most twice the model gap + 300 bp."""
    out: list[Hit] = [h for h in hits if h.region not in tiled]
    hits = [h for h in hits if h.region in tiled]
    for (_, strand), hs in itertools.groupby(sorted(hits, key=lambda h: (h.region, h.strand, h.start)),
                                             key=lambda h: (h.region, h.strand)):
        chains: list[list[Hit]] = []
        for h in hs:
            for c in chains:
                p = c[-1]
                if strand == "+":
                    ahead, model_gap = h.hmm_to > p.hmm_to and h.hmm_from >= p.hmm_from, h.hmm_from - p.hmm_to
                else:  # the model runs backwards along the record
                    ahead, model_gap = h.hmm_from < p.hmm_from and h.hmm_to <= p.hmm_to, p.hmm_from - h.hmm_to
                if ahead and h.end > p.end and h.start - p.end <= 2 * max(model_gap, 0) + 300:
                    c.append(h)
                    break
            else:
                chains.append([h])
        for c in chains:
            if len(c) == 1:
                out.append(c[0])
                continue
            covered = set().union(*(range(h.start, h.end + 1) for h in c))
            # tiles overlap: count each residue's score once (summed score scaled to the covered length)
            score = sum(h.score for h in c) * len(covered) / sum(h.end - h.start + 1 for h in c)
            out.append(Hit(c[0].region, min(h.start for h in c), max(h.end for h in c),
                           min(h.hmm_from for h in c), max(h.hmm_to for h in c), c[0].hmm_len, round(score, 1), strand))
    return out


def verify(inputs: list[str], g2t_dir: str, out_dir: str, models: list[Path] | None = None,
           all_records: bool = False, cpus: int = 0, min_bp: int = 150, min_cover: float = 0.5,
           cache: Path | None = None, cascade: bool = True) -> pd.DataFrame:
    """``cache``: SQLite file of hits per sequence and model (None: no cache); ``cascade``: run the
    models of g2t's region first, the rest only for records they do not explain."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    status = read_status(g2t_dir)
    keep = select(status, all_records)
    gtype = dict(zip(status["acc_base"], status["gene_type"]))
    stat = dict(zip(status["acc_base"], status["status"]))
    recs = {a: (o, s) for a, o, s in read_sequences(inputs, keep)}
    logger.info("scanning %d of %d records", len(recs), len(status))
    whole = load_models(models or default_models())
    ms = tile_models(whole)
    hc = HitCache(cache)
    seqs = {a: s for a, (_, s) in recs.items()}
    every = {m.region for m in ms}
    want = {a: (stage1_regions(g2t_types(gtype.get(a.split(".")[0]))) if cascade else None) or every
            for a in seqs}
    hits = scan(seqs, ms, cpus=cpus, want=want, cache=hc) if seqs else {}
    rest = {a: every - w for a, w in want.items() if w != every and not explained(hits.get(a, []), len(seqs[a]))}
    if rest:
        logger.info("%d of %d records not explained by the models of their g2t type: other models",
                    len(rest), len(seqs))
        for a, h in scan({a: seqs[a] for a in rest}, ms, cpus=cpus, want=rest, cache=hc).items():
            hits.setdefault(a, []).extend(h)
    # tile ends that stop short of the neighbouring region / record end: the full model, there only
    tiled = {m.region for m in ms if m.full}
    ends = {f"{a}\t{k}\t{x}": seqs[a][x - 1: y] for a, hs in hits.items()
            for k, x, y in ragged(hs, len(seqs[a]), tiled)}
    if ends:
        logger.info("%d region ends the tiles did not reach: full model on those stretches", len(ends))
        got = scan(ends, whole, cpus=cpus, want={e: {e.split("\t")[1]} for e in ends}, cache=hc)
        for e, hs in got.items():
            a, k, first = e.split("\t")
            shift = int(first) - 1
            hits[a] = extend(hits[a], k, [Hit(h.region, h.start + shift, h.end + shift, h.hmm_from, h.hmm_to,
                                              h.hmm_len, h.score, h.strand) for h in hs if h.region == k])
    # records without any hit: once more at E <= 1e-3 (weak partial 5.8S of ITS records)
    again = {a: s for a, s in seqs.items() if a not in hits}
    if again:
        logger.info("rescanning %d records without a hit at E <= %g", len(again), RESCUE_EVALUE)
        hits.update(scan(again, ms, cpus=cpus, evalue=RESCUE_EVALUE, want={a: want[a] for a in again}, cache=hc))
    rows = []
    fasta: dict[str, list[str]] = {}
    for a, (org, s) in recs.items():
        r = regions(hits.get(a, []), len(s))
        sup = supported_types(r, min_bp, min_cover)
        g = g2t_types(gtype.get(a.split(".")[0]))
        v = verdict(g, sup, bool(r.spans))
        rows.append({"accession": a, "organism": org, "length": len(s), "g2t_status": stat.get(a.split(".")[0]),
                     "g2t_type": ",".join(g), "sequence_type": ",".join(sup), "verdict": v,
                     "regions": _fmt(r), "strand": r.strand if r.spans else ""})
        for k, (x, y, _) in r.spans.items():
            seg = s[x - 1: y]
            if r.strand == "-":
                from Bio.Seq import Seq

                seg = str(Seq(seg).reverse_complement())
            fasta.setdefault(k, []).append(f">{a} {org} {k}:{x}-{y}\n{seg}\n")
    t = pd.DataFrame(rows, columns=["accession", "organism", "length", "g2t_status", "g2t_type", "sequence_type",
                                    "verdict", "regions", "strand"])
    t.to_csv(out / "verify.csv", index=False)
    rd = out / "regions"
    rd.mkdir(exist_ok=True)
    for k, lines in fasta.items():
        (rd / f"{k}.fa").write_text("".join(lines), encoding="utf-8")
    return t


def describe(t: pd.DataFrame) -> list[str]:
    lines = [f"{len(t)} records scanned: " + ", ".join(f"{k} {v}" for k, v in t["verdict"].value_counts().items())]
    for v in ("different", "new", "no_hit"):
        d = t[t["verdict"] == v]
        if len(d):
            pairs = (d["g2t_type"].replace("", "-") + " -> " + d["sequence_type"].replace("", "-")).value_counts()
            lines.append(f"{v}: " + ", ".join(f"{k} ({n})" for k, n in pairs.head(6).items()))
    return lines


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="g2t-verify", description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog="examples:\n  g2t-verify -i gb/markers --g2t out -o verify\n"
                                       "  g2t-verify -i gb/markers --g2t out -o verify --models my_rrna.hmm")
    p.add_argument("-i", "--input", nargs="+", required=True, help="GenBank files or folders (with sequences)")
    p.add_argument("--g2t", required=True, metavar="DIR", help="output folder of the g2t run on these records")
    p.add_argument("-o", "--output", required=True, metavar="DIR")
    p.add_argument("--models", nargs="+", metavar="HMM", help="your own nucleotide HMM file(s) instead of barrnap 0.9")
    p.add_argument("--all", action="store_true", help="scan every record (default: rRNA, unmatched, conflicted)")
    p.add_argument("--min-bp", type=int, default=150, help="region length that counts as a gene (default 150)")
    p.add_argument("--min-cover", type=float, default=0.5, help="or this share of the model (default 0.5)")
    p.add_argument("--cpu", type=int, default=0, help="threads (default: all cores)")
    p.add_argument("--no-cache", action="store_true",
                   help=f"do not reuse or store hits (default cache: {hit_cache_path()})")
    p.add_argument("--all-models", action="store_true",
                   help="run every model on every record (default: g2t's region first, others if needed)")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    import time

    t0 = time.monotonic()
    t = verify(a.input, a.g2t, a.output, [Path(m) for m in a.models] if a.models else None, a.all, a.cpu,
               a.min_bp, a.min_cover, cache=None if a.no_cache else hit_cache_path(), cascade=not a.all_models)
    print(f"time: {time.monotonic() - t0:.1f} s")
    for line in describe(t):
        print(line)
    print(f"output: {a.output}/verify.csv, {a.output}/regions/*.fa")
    return 0


if __name__ == "__main__":
    sys.exit(main())
