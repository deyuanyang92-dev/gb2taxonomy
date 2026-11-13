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
from 50 bp. Records without any hit are scanned a second time without acceleration filters.
Overlapping hits of different models (mitochondrial 16S also matches nuclear LSU) compete: the
higher score wins.

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


def scan(seqs: dict[str, str], model_files: list[Path], cpus: int = 0, sensitive: bool = False
         ) -> dict[str, list[Hit]]:
    """nhmmer of every model against the sequences (pyhmmer); hits with E <= 1e-5.

    ``sensitive``: no acceleration filters (nhmmer --max), for short fragments the filters drop."""
    try:
        import pyhmmer
        from pyhmmer.easel import DigitalSequenceBlock, TextSequence
        from pyhmmer.plan7 import HMMFile
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise SystemExit("g2t-verify needs pyhmmer: pip install pyhmmer  (or pip install -e '.[verify]')") from e
    hmms = []
    for f in model_files:
        with HMMFile(str(f)) as hf:
            hmms += [h for h in hf if model_region(_s(h.name))]
    if not hmms:
        raise SystemExit("no 12S/16S/18S/5.8S/28S models in " + ", ".join(map(str, model_files)))
    alphabet = hmms[0].alphabet
    names = list(seqs)
    block = DigitalSequenceBlock(alphabet, [
        TextSequence(name=n.encode(), sequence=_clean(seqs[n])).digitize(alphabet) for n in names])
    out: dict[str, list[Hit]] = {}
    opts: dict[str, Any] = dict(F1=1.0, F2=1.0, F3=1.0, bias_filter=False) if sensitive else {}
    for hits in pyhmmer.hmmer.nhmmer(hmms, block, cpus=cpus, E=EVALUE, **opts):
        q = hits.query
        region = model_region(_s(q.name))
        for h in hits:
            if not h.included:
                continue
            for d in h.domains:
                if not d.included:
                    continue
                a, b = sorted((d.env_from, d.env_to))
                out.setdefault(_s(h.name), []).append(Hit(
                    region or "", a, b, d.alignment.hmm_from, d.alignment.hmm_to, q.M, d.score,
                    "+" if d.env_from <= d.env_to else "-"))
    return out


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


def read_sequences(paths: Iterable[str | Path], keep: set[str] | None = None) -> Iterator[tuple[str, str, str]]:
    """(accession.version, organism, sequence) of GenBank records; records without sequence are skipped."""
    import warnings

    from Bio import BiopythonParserWarning, SeqIO
    from Bio.Seq import UndefinedSequenceError

    from g2t._extract.files import open_genbank
    from g2t.mitostats import _gb_files

    for f in _gb_files(paths):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open_genbank(str(f)) as fh:
                for rec in SeqIO.parse(fh, "genbank"):
                    if keep is not None and rec.id.split(".")[0] not in keep:
                        continue
                    try:
                        s = str(rec.seq)
                    except UndefinedSequenceError:  # CONTIG / master records
                        continue
                    if s:
                        yield rec.id, rec.annotations.get("organism", ""), s


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


def verify(inputs: list[str], g2t_dir: str, out_dir: str, models: list[Path] | None = None,
           all_records: bool = False, cpus: int = 0, min_bp: int = 150, min_cover: float = 0.5) -> pd.DataFrame:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    status = read_status(g2t_dir)
    keep = select(status, all_records)
    gtype = dict(zip(status["acc_base"], status["gene_type"]))
    stat = dict(zip(status["acc_base"], status["status"]))
    recs = {a: (o, s) for a, o, s in read_sequences(inputs, keep)}
    logger.info("scanning %d of %d records", len(recs), len(status))
    model_files = models or default_models()
    hits = scan({a: s for a, (_, s) in recs.items()}, model_files, cpus=cpus) if recs else {}
    # second pass without acceleration filters for records without any hit (short 28S fragments)
    again = {a: s for a, (_, s) in recs.items() if a not in hits}
    if again:
        logger.info("rescanning %d records without a hit (sensitive mode)", len(again))
        hits.update(scan(again, model_files, cpus=cpus, sensitive=True))
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
    p.add_argument("--cpu", type=int, default=0, help="threads (default: all)")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    t = verify(a.input, a.g2t, a.output, [Path(m) for m in a.models] if a.models else None, a.all, a.cpu,
               a.min_bp, a.min_cover)
    for line in describe(t):
        print(line)
    print(f"output: {a.output}/verify.csv, {a.output}/regions/*.fa")
    return 0


if __name__ == "__main__":
    sys.exit(main())
