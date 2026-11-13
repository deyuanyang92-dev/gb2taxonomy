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
    same_species_far  two records with the same species name are > 15% apart in a fast-evolving
                      marker: one is misidentified or contaminated, or the name covers cryptic
                      species. The note gives each record's distance to its own genus: the record far
                      from its genus is the suspect. E.g. HM746718 "Ophelia neglecta" 16S is 24.3%
                      from the 16S of the O. neglecta mitogenome OZ260669 and 23% from every other
                      Ophelia, but 16% from Ophelina. Genes of mitogenomes are included (12S / 16S
                      found by the g2t-verify HMMs when the record has no annotation).
    remote            no sequence of this marker in the data set within a k-mer distance of 0.35
                      (roughly > 30% different): wrong gene, wrong organism, or the only sample of
                      a distant lineage. Data sets with fewer than 10 sequences are not judged.

A flag is a prompt to look, not a verdict: a genus may be para- or polyphyletic, a name may be
outdated, and a sister genus can share haplotypes. Check flagged records against the paper,
the voucher and a tree.

Method (alignment-free, no external programs): canonical k-mers (both strands, so reversed
records and introns do not break it). Candidate neighbours by cosine similarity of 5-mer counts;
distances from 11-mers shared with the shorter sequence (1 - C^(1/k), so a fragment is not penalised
against a complete gene), and the Mash distance on 8-mers for the remote test.
An earlier alignment-based prototype (MAFFT + p-distance) put AY364856 Ophelia
rathkei and LC670764 Travisia at 49-51% from everything; both are ~21-25% from their own genus
by k-mers: MAFFT misaligned a reversed / intron-containing COI. See docs/evaluation/screen.md.

CLI:
    g2t-screen -i gb/markers --g2t out              # writes out/organized_genes/screen.csv
"""

from __future__ import annotations

import argparse
import contextlib
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
SAME_FAR = 0.15  # same_species_far: two records of one named species > 15% apart (fast markers)
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


def pdist(a: np.ndarray, b: np.ndarray, k: int) -> float:
    """p-distance estimate from the shared k-mers of the shorter sequence: 1 - C^(1/k), C = |a & b| /
    min(|a|, |b|). Unlike Jaccard / Mash it does not grow when a fragment is compared with a longer
    sequence of the same gene (a 544 bp 16S against a mitogenome's 1,231 bp 16S)."""
    m = min(a.size, b.size)
    if m == 0:
        return 1.0
    c = np.intersect1d(a, b, assume_unique=True).size / m
    return 1.0 if c == 0 else max(0.0, 1.0 - c ** (1.0 / k))


def species_of(organism: str) -> str:
    """'Genus epithet' for a named species; '' for sp. / cf. / aff. / complex / environmental names."""
    w = (organism or "").split()
    if len(w) < 2 or not w[0][:1].isupper() or not re.fullmatch(r"[a-z][a-z-]+", w[1]) \
            or re.search(r"\b(sp|spp|cf|aff|nr|complex|environmental|sample|uncultured)\b", organism):
        return ""
    return f"{w[0]} {w[1]}"


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


def screen_marker(ids: list[str], seqs: list[str], genera: list[str], fast: bool,
                  species: list[str] | None = None) -> list[dict]:
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
    species = species or [""] * n
    by_species: dict[str, list[int]] = {}
    for i, sp in enumerate(species):
        if sp:
            by_species.setdefault(sp, []).append(i)
    rows = []
    for start in range(0, n, 2000):  # cosine in row blocks: memory stays ~ 2000 x n
        c = x[start : start + 2000] @ x.T
        for r, i in enumerate(range(start, min(n, start + 2000))):
            c[r, i] = -1
            cand = np.argsort(-c[r])[: min(CANDIDATES, n - 1)]
            d11 = sorted((pdist(k11[i], k11[j], 11), int(j)) for j in cand)
            near_d, near_j = d11[0]
            remote_d = min(mash(k8[i], k8[j], 8) for j in cand)
            g, gj = genera[i], genera[near_j]
            flag, note = "", ""
            mates = [j for j in by_genus.get(g, []) if j != i] if g else []
            if len(mates) > 5:  # the 5 most similar of its genus by 5-mer profile
                mates = [mates[t] for t in np.argsort(-c[r, mates])[:5]]
            own = [pdist(k11[i], k11[j], 11) for j in mates]
            same = [(pdist(k11[i], k11[j], 11), j) for j in by_species.get(species[i], []) if j != i]
            same_d, same_j = min(same) if same else (None, -1)
            own_d = min(own) if own else None
            if fast and g and gj and gj != g and near_d <= NEAR and (own_d is None or own_d > near_d + MARGIN):
                flag = "near_other_genus"
                note = f"{near_d:.1%} from {gj}" + (
                    f"; nearest {g} {own_d:.1%}" if own_d is not None else f"; no other {g} sequence")
            elif fast and same_d is not None and same_d > SAME_FAR:
                flag = "same_species_far"
                own_txt = f"nearest {g} {own_d:.1%}" if own_d is not None else f"no other {g}"
                note = (f"{same_d:.1%} from {ids[same_j]}, also named {species[i]}; {own_txt}; "
                        f"nearest overall {gj or '?'} {near_d:.1%}")
            elif n >= MIN_SET and remote_d > REMOTE:
                flag, note = "remote", f"nearest sequence {remote_d:.2f} (8-mer distance)"
            if flag:
                rows.append({"accession": ids[i], "genus": g, "flag": flag, "note": note,
                             "nearest": ids[near_j], "nearest_genus": gj, "distance": round(near_d, 4),
                             "own_genus_distance": None if own_d is None else round(own_d, 4)})
    return rows


_MT_GENES = [("coi", r"\bco(x)?1\b|\bcoi\b|cytochrome c oxidase subunit (1|i)\b"),
             ("cox2", r"\bco(x)?2\b|\bcoii\b|cytochrome c oxidase subunit (2|ii)\b"),
             ("cox3", r"\bco(x)?3\b|\bcoiii\b|cytochrome c oxidase subunit (3|iii)\b"),
             ("cob", r"\bcob\b|\bcytb\b|cytochrome b\b"),
             ("nd4l", r"\bnad?4l\b|nadh dehydrogenase subunit 4l"),
             *[(f"nd{i}", rf"\bnad?{i}\b|nadh dehydrogenase subunit {i}\b") for i in (1, 2, 3, 4, 5, 6)],
             ("atp6", r"\batp(ase)? ?6\b|atp synthase f0 subunit 6"),
             ("atp8", r"\batp(ase)? ?8\b|atp synthase f0 subunit 8"),
             ("16s", r"\b16s\b|\brrnl\b|\bl-rrna\b|large subunit ribosomal"),
             ("12s", r"\b12s\b|\brrns\b|\bs-rrna\b|small subunit ribosomal")]


def mitogenome_genes(paths: Iterable[str], keep: set[str]) -> dict[str, dict[str, str]]:
    """Genes of mitogenome records: {accession.version: {gene: sequence}}. Protein genes from the
    annotation (first CDS of each gene); 12S / 16S from the annotation or, where it has none, by the
    g2t-verify profile HMMs. A 16S fragment can then be compared with the 16S of a mitogenome of the
    same species (HM746718 "Ophelia neglecta" vs the O. neglecta mitogenome OZ260669: 24% apart)."""
    import warnings

    from Bio import BiopythonParserWarning, SeqIO

    from g2t._extract.files import open_genbank
    from g2t.mitostats import _gb_files

    out: dict[str, dict[str, str]] = {}
    raw: dict[str, str] = {}  # records whose 12S / 16S the annotation does not give
    pats = [(g, re.compile(p, re.I)) for g, p in _MT_GENES]
    for f in _gb_files(paths):
        with warnings.catch_warnings(), open_genbank(str(f)) as fh:
            warnings.simplefilter("ignore", BiopythonParserWarning)
            for rec in SeqIO.parse(fh, "genbank"):
                if rec.id.split(".")[0] not in keep:
                    continue
                genes: dict[str, str] = {}
                for ft in rec.features:
                    if ft.type not in ("CDS", "rRNA"):
                        continue
                    label = " ".join(ft.qualifiers.get("gene", []) + ft.qualifiers.get("product", []))
                    for g, p in pats:
                        if g not in genes and p.search(label) and (g in ("12s", "16s")) == (ft.type == "rRNA"):
                            with contextlib.suppress(Exception):  # undefined sequence etc.
                                genes[g] = str(ft.extract(rec.seq)).upper()
                            break
                if genes:
                    out[rec.id] = genes
                if not {"12s", "16s"} <= set(genes):
                    with contextlib.suppress(Exception):
                        raw[rec.id] = str(rec.seq).upper()
    # 12S / 16S by profile HMM where the annotation has none or unnamed rRNA features (OZ260669 has no
    # features at all): the same models as g2t-verify; skipped without pyhmmer / models
    if raw:
        with contextlib.suppress(Exception, SystemExit):
            from g2t import verify as vf

            models = [m for m in vf.tile_models(vf.load_models(vf.default_models())) if m.region in ("12s", "16s")]
            hits = vf.scan(raw, models, cache=vf.HitCache(vf.hit_cache_path()))
            for acc, hs in hits.items():
                r = vf.regions(hs, len(raw[acc]))
                for g in ("12s", "16s"):
                    if g in r.spans and g not in out.get(acc, {}):
                        x, y, _ = r.spans[g]
                        seg = raw[acc][x - 1 : y]
                        if r.strand == "-":
                            seg = seg.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                        out.setdefault(acc, {})[g] = seg
    return out


def screen(inputs: Iterable[str], g2t_dir: str | Path) -> pd.DataFrame:
    """Screen every marker of a g2t run (labeled_genes/assigned_genes_types_all.csv); genes of
    mitogenome records are taken from their annotation and screened with the marker records."""
    from g2t.verify import read_sequences

    inputs = list(inputs)
    d = Path(g2t_dir)
    a = pd.read_csv(d / "labeled_genes" / "assigned_genes_types_all.csv", dtype=str, keep_default_na=False,
                    low_memory=False)
    oname = "Organism" if "Organism" in a else "organism"
    org = dict(zip(a["ACCESSION"], a[oname]))
    gen = dict(zip(a["ACCESSION"], a["Genus"] if "Genus" in a else [""] * len(a)))
    mt = set(a.loc[a["gene_type"] == "mtgenome", "ACCESSION"])
    a = a[~a["gene_type"].str.contains(",") & (a["gene_type"] != "") & ~a["gene_type"].isin(SKIP)]
    seqs = {acc.split(".")[0]: (acc, s) for acc, _, s in read_sequences(inputs, set(a["ACCESSION"]))}
    by_gene: dict[str, list[tuple[str, str]]] = {}  # gene -> [(accession.version, sequence)]
    for acc, gene in zip(a["ACCESSION"], a["gene_type"]):
        if acc in seqs:
            by_gene.setdefault(gene, []).append(seqs[acc])
    for acc, genes in (mitogenome_genes(inputs, mt) if mt else {}).items():
        for gene, s in genes.items():
            by_gene.setdefault(gene, []).append((acc, s))
    cols = ["gene", "accession", "organism", "genus", "flag", "note", "nearest", "nearest_organism", "nearest_genus",
            "distance", "own_genus_distance"]
    out = []
    for gene, items in sorted(by_gene.items()):
        if len(items) < 2:
            continue
        ids = [x for x, _ in items]
        bases = [x.split(".")[0] for x in ids]
        rows = screen_marker(ids, [s for _, s in items], [genus_of(org.get(b, ""), gen.get(b, "")) for b in bases],
                             gene in FAST, [species_of(org.get(b, "")) for b in bases])
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
    genes = [c for c in m.columns if c in set(flags["gene"]) | {"mtgenome"}] if len(flags) else []

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
