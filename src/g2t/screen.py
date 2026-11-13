"""g2t-screen: flag possible contaminated or misidentified sequences, marker by marker.

GenBank holds sequences whose organism is wrong: contamination (another animal's DNA, often a
lab-mate's), misidentified specimens, or sample mix-ups. In a phylogeny they show up as a species
nested inside another genus. Checking this by hand (align, compute distances, build a tree, look)
is slow; g2t-screen does the distance part for every marker in seconds and flags:

    near_other_genus  the closest sequence is of another genus and almost identical (<= 2% k-mer
                      distance), while the closest sequence of its own genus (if any) is clearly
                      farther. E.g. Opheliidae COI: OQ417162 "Polyophthalmus pictus" is 0.4% from
                      an Armandia sequence. Only for fast-evolving markers (COI, 16S, 12S, CYTB,
                      ND genes, ITS ...): 18S / 28S / H3 can be identical across genera.
    remote            no sequence of this marker in the data set within a k-mer distance of 0.35
                      (roughly > 30% different): wrong gene, wrong organism, or the only sample of
                      a distant lineage. Data sets with fewer than 10 sequences are not judged.

A flag is a prompt to look, not a verdict: a genus may be para- or polyphyletic, a name may be
outdated, and a sister genus can share haplotypes. Check flagged records against the paper,
the voucher and a tree.

Method (alignment-free, no external programs): canonical k-mers (both strands, so reversed
records and introns do not break it). Candidate neighbours by cosine similarity of 5-mer counts,
then the Mash distance -ln(2J/(1+J))/k on 11-mers (accurate for small distances) and 8-mers (for
the remote test). An earlier alignment-based prototype (MAFFT + p-distance) put AY364856 Ophelia
rathkei and LC670764 Travisia at 49-51% from everything; both are ~21-25% from their own genus
by k-mers: MAFFT misaligned a reversed / intron-containing COI. See docs/evaluation/screen.md.

CLI:
    g2t-screen -i gb/markers --g2t out              # writes out/organized_genes/screen.csv
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger("g2t.screen")

FAST = {"coi", "cox2", "cox3", "16s", "12s", "cob", "nd1", "nd2", "nd3", "nd4", "nd4l", "nd5", "nd6", "atp6", "atp8",
        "its1-its2", "matk", "rbcl", "trnl-trnf", "psba-trnh"}
SKIP = {"mtgenome", "plastome", "18-28s"}
NEAR = 0.02  # k-mer (11-mer Mash) distance to another genus
MARGIN = 0.02  # own genus must be at least this much farther
REMOTE = 0.35  # 8-mer Mash distance to the nearest sequence
MIN_SET = 10
CANDIDATES = 20


_LUT = np.full(256, 4, np.uint8)
for _i, _c in enumerate(b"ACGT"):
    _LUT[_c] = _i


def _windows(s: str, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Forward and reverse-complement 2-bit codes of every ACGT-only k-mer window."""
    c = _LUT[np.frombuffer(s.encode("ascii", "replace"), np.uint8)]
    n = len(c) - k + 1
    if n <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, np.uint64)
    bad = np.concatenate([[0], np.cumsum(c == 4)])
    ok = (bad[k:] - bad[:-k]) == 0
    c64 = np.where(c == 4, 0, c).astype(np.uint64)
    fwd = np.zeros(n, np.uint64)
    rev = np.zeros(n, np.uint64)
    for j in range(k):
        fwd = (fwd << np.uint64(2)) | c64[j : j + n]
        rev = rev | ((np.uint64(3) - c64[j : j + n]) << np.uint64(2 * j))
    return fwd[ok], rev[ok]


def kmers(s: str, k: int) -> np.ndarray:
    """Sorted unique canonical k-mers (both strands; windows with N etc. skipped), 2-bit coded."""
    fwd, rev = _windows(s, k)
    return np.unique(np.minimum(fwd, rev))


def mash(a: np.ndarray, b: np.ndarray, k: int) -> float:
    """Mash distance -ln(2J/(1+J))/k of two k-mer sets (0 for identical sets)."""
    inter = np.intersect1d(a, b, assume_unique=True).size
    u = a.size + b.size - inter
    j = inter / u if u else 0.0
    return 1.0 if j == 0 else max(0.0, float(-np.log(2 * j / (1 + j)) / k))


def profile(seqs: list[str]) -> np.ndarray:
    """Unit-length 5-mer count vectors (both strands)."""
    x = np.zeros((len(seqs), 1024), np.float32)
    for n, s in enumerate(seqs):
        fwd, rev = _windows(s, 5)
        x[n] = np.bincount(fwd.astype(np.int64), minlength=1024) + np.bincount(rev.astype(np.int64), minlength=1024)
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.where(norm == 0, 1, norm)


def genus_of(organism: str, genus: str = "") -> str:
    """The genus, or '' when unresolved ("Opheliidae sp.", "uncultured ...", empty)."""
    g = (genus or "").strip() or (organism or "").split(" ")[0]
    if not g or not g[0].isupper() or g.lower() in ("uncultured", "unidentified", "environmental") \
            or re.search(r"(idae|inae|ini|ales|aceae|formes|ida|ea)$", g):
        return ""
    return g


def screen_marker(ids: list[str], seqs: list[str], genera: list[str], fast: bool) -> list[dict]:
    """Flags for one marker's sequences (same order); one row per flagged sequence."""
    n = len(ids)
    if n < 2:
        return []
    x = profile(seqs)
    k11 = [kmers(s, 11) for s in seqs]
    k8 = [kmers(s, 8) for s in seqs]
    by_genus: dict[str, list[int]] = {}
    for i, g in enumerate(genera):
        if g:
            by_genus.setdefault(g, []).append(i)
    rows = []
    for start in range(0, n, 2000):  # cosine in row blocks: memory stays ~ 2000 x n
        c = x[start : start + 2000] @ x.T
        for r, i in enumerate(range(start, min(n, start + 2000))):
            c[r, i] = -1
            cand = np.argsort(-c[r])[: min(CANDIDATES, n - 1)]
            d11 = sorted((mash(k11[i], k11[j], 11), int(j)) for j in cand)
            near_d, near_j = d11[0]
            remote_d = min(mash(k8[i], k8[j], 8) for j in cand)
            g, gj = genera[i], genera[near_j]
            flag, note = "", ""
            mates = [j for j in by_genus.get(g, []) if j != i] if g else []
            if len(mates) > 5:  # the 5 most similar of its genus by 5-mer profile
                mates = [mates[t] for t in np.argsort(-c[r, mates])[:5]]
            own = [mash(k11[i], k11[j], 11) for j in mates]
            own_d = min(own) if own else None
            if fast and g and gj and gj != g and near_d <= NEAR and (own_d is None or own_d > near_d + MARGIN):
                flag = "near_other_genus"
                note = f"{near_d:.1%} from {gj}" + (
                    f"; nearest {g} {own_d:.1%}" if own_d is not None else f"; no other {g} sequence")
            elif n >= MIN_SET and remote_d > REMOTE:
                flag, note = "remote", f"nearest sequence {remote_d:.2f} (8-mer distance)"
            if flag:
                rows.append({"accession": ids[i], "genus": g, "flag": flag, "note": note,
                             "nearest": ids[near_j], "nearest_genus": gj, "distance": round(near_d, 4),
                             "own_genus_distance": None if own_d is None else round(own_d, 4)})
    return rows


def screen(inputs: Iterable[str], g2t_dir: str | Path) -> pd.DataFrame:
    """Screen every marker of a g2t run (labeled_genes/assigned_genes_types_all.csv)."""
    from g2t.verify import read_sequences

    d = Path(g2t_dir)
    a = pd.read_csv(d / "labeled_genes" / "assigned_genes_types_all.csv", dtype=str, keep_default_na=False,
                    low_memory=False)
    a = a[~a["gene_type"].str.contains(",") & (a["gene_type"] != "") & ~a["gene_type"].isin(SKIP)]
    keep = set(a["ACCESSION"])
    seq = {acc.split(".")[0]: (acc, s) for acc, _, s in read_sequences(inputs, keep)}
    cols = ["gene", "accession", "organism", "genus", "flag", "note", "nearest", "nearest_organism", "nearest_genus",
            "distance", "own_genus_distance"]
    out = []
    org = dict(zip(a["ACCESSION"], a.get("Organism", a.get("organism", ""))))
    for gene, g in a.groupby("gene_type"):
        g = g[g["ACCESSION"].isin(seq)]
        if len(g) < 2:
            continue
        ids = [seq[x][0] for x in g["ACCESSION"]]
        genera = [genus_of(o, ge) for o, ge in zip(g.get("Organism", g.get("organism")), g.get("Genus", [""] * len(g)))]
        rows = screen_marker(ids, [seq[x][1] for x in g["ACCESSION"]], genera, gene in FAST)
        for r in rows:
            r.update(gene=gene, organism=org.get(r["accession"].split(".")[0], ""),
                     nearest_organism=org.get(r["nearest"].split(".")[0], ""))
        out += rows
        logger.info("%s: %d sequences, %d flagged", gene, len(ids), len(rows))
    return pd.DataFrame(out, columns=cols)


def mark_matrix(matrix_csv: str | Path, flags: pd.DataFrame) -> int:
    """Add a ``screen_flags`` column to the specimen x gene matrix; returns rows flagged."""
    m = pd.read_csv(matrix_csv, dtype=str, keep_default_na=False, low_memory=False)
    text = {r.accession: f"{r.accession} {r.gene}: {r.flag} ({r.note})" for r in flags.itertuples()}
    genes = [c for c in m.columns if c in set(flags["gene"])] if len(flags) else []

    def row_flags(r: pd.Series) -> str:
        accs = [x.strip() for c in genes for x in re.split(r"[;,|]", str(r[c])) if x.strip()]
        return " | ".join(text[x] for x in accs if x in text)

    m["screen_flags"] = m.apply(row_flags, axis=1) if len(m) else []
    m.to_csv(matrix_csv, index=False)
    return int((m["screen_flags"] != "").sum())


def pipeline_step(inputs: list[str], out_root: str | Path) -> list[str]:
    """Step 5 of ``g2t``: screen every marker, write organized_genes/screen.csv, mark the matrix."""
    root = Path(out_root)
    t = screen(inputs, root)
    t.to_csv(root / "organized_genes" / "screen.csv", index=False)
    matrix = root / "organized_genes" / "organized_species_voucher.csv"
    n = mark_matrix(matrix, t) if matrix.exists() else 0
    counts = t["flag"].value_counts().to_dict() if len(t) else {}
    return [f"Step 5: contamination / misidentification screen: {len(t)} sequences flagged"
            + (" (" + ", ".join(f"{k} {v}" for k, v in counts.items()) + ")" if counts else "")
            + f" in {n} specimens; organized_genes/screen.csv, column screen_flags"]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="g2t-screen", description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", nargs="+", required=True, help="GenBank files or folders (with sequences)")
    p.add_argument("--g2t", required=True, metavar="DIR", help="output folder of the g2t run on these records")
    p.add_argument("-o", "--output", metavar="CSV", help="default: <g2t>/organized_genes/screen.csv")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    t = screen(a.input, a.g2t)
    out = Path(a.output) if a.output else Path(a.g2t) / "organized_genes" / "screen.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(out, index=False)
    print(f"{len(t)} sequences flagged; {out}")
    if len(t):
        print(t[["gene", "accession", "organism", "flag", "note"]].head(30).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
