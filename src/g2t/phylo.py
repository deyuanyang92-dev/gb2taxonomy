"""Phylogenetic screen: does each sequence sit where its name says, in a gene tree with relatives?

Design (cf. OrthoFinder 3, Emms et al. 2026, Nat. Methods 23:1327: gather over-inclusively, then let a
gene tree decide membership; FAMSA alignment, gappy columns trimmed, FastTree):

1. Per marker, the study taxon's sequences (all of them, or one per species plus every sequence the
   distance screen flagged when there are more than ``cap``) and background sequences of related
   families (one per species, at most ``bg_per_family`` per family) - the families where a
   contaminant or a misplaced name most likely belongs.
2. Sequences are put on one strand (k-mers against the longest sequences; FAMSA does not reverse),
   aligned with FAMSA (MAFFT if FAMSA is missing), columns with >= 90% gaps dropped, tree by FastTree
   (GTR+G; -fastest above 1,000 sequences) or IQ-TREE (ModelFinder, 1,000 ultrafast bootstraps).
3. Unrooted test: for each study sequence, the smallest side of a well-supported split (support >=
   ``support``) that holds it and >= 2 other well-covered sequences = its neighbours. No root is
   needed, so the result does not depend on which background family is the outgroup.
4. Record or taxon? A sequence is flagged when its neighbours are another family (all of them, one
   or two families) or other genera of its family, while other sequences of its own genus sit
   elsewhere: the record does not match its name (contamination, mix-up, misidentification). When
   its congeners sit with it, or more than half of a genus is flagged, the genus as a whole is placed
   there: a taxonomic question (non-monophyletic genus, outdated name), reported separately. A
   sequence that is the only one of its genus cannot be told apart and says so.
5. Sequences covering < ``min_cover`` of the core region of the alignment (columns where most
   sequences have data; a record of another region of the gene, e.g. a different 28S D-region)
   and neighbours mixing the own family with others are not judged.
6. Distance cross-check: a flag stands only if the sequence is also closer (shared k-mers) to the
   group it sits in than to its nearest congener; a close congener vetoes a placement artefact.
7. Cross-marker check per specimen (g2t matrix): a misplaced COI whose 18S / 28S sit correctly is a
   bad sequence (contamination or a mix-up of that sequence); all markers of a specimen misplaced
   means the specimen is misidentified.

Opheliidae (2026-11-13): HM746718 "Ophelia neglecta" 16S sits in Notomastus (Capitellidae, 0.99);
PX620130 "O. translucens" 16S in Hyboscolex (Scalibregmatidae, 0.99; 0.2% from a Hyboscolex);
OQ417162 "Polyophthalmus pictus" COI in Armandia (1.00). See docs/evaluation/phylo-screen.md.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from g2t import screen as sc

logger = logging.getLogger("g2t.phylo")

NEAR_FAMILY = 0.05  # a lone sequence this close to another family is reported as a contaminant
MARKERS = ("coi", "16s", "12s", "18s", "28s", "cob", "h3", "its1-its2")
_COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


# ------------------------------------------------------------------ tools


def find_tool(*names: str) -> str | None:
    """First executable found on PATH or in the conda environments of this machine."""
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    base = Path(os.environ.get("CONDA_PREFIX_1") or os.environ.get("CONDA_PREFIX") or Path.home() / "miniconda3")
    for root in (base, base.parent):
        for d in sorted(root.glob("envs/*/bin")) + [root / "bin"]:
            for n in names:
                if (d / n).exists():
                    return str(d / n)
    return None


@dataclass
class Tools:
    aligner: str
    aligner_kind: str  # famsa | mafft
    tree: str
    tree_kind: str  # fasttree | iqtree

    @classmethod
    def find(cls, iqtree: bool = False) -> Tools:
        famsa, mafft = find_tool("famsa"), find_tool("mafft")
        ft, iq = find_tool("FastTree", "fasttree", "FastTreeMP"), find_tool("iqtree3", "iqtree2", "iqtree")
        if not (famsa or mafft):
            raise SystemExit("g2t-screen --tree needs FAMSA or MAFFT (conda install -c bioconda famsa fasttree)")
        tree, kind = (iq, "iqtree") if iqtree and iq else (ft, "fasttree") if ft else (iq, "iqtree")
        if not tree:
            raise SystemExit("g2t-screen --tree needs FastTree or IQ-TREE (conda install -c bioconda fasttree)")
        return cls(famsa or mafft or "", "famsa" if famsa else "mafft", tree, kind)


# ------------------------------------------------------------------ alignment and tree


def orient(seqs: dict[str, str], n_ref: int = 5) -> tuple[dict[str, str], list[str]]:
    """Put every sequence on the strand of the ``n_ref`` longest ones (shared 11-mers)."""
    refs = sorted(seqs, key=lambda k: len(seqs[k]), reverse=True)[:n_ref + 1]
    fwd = {k: np.unique(sc._windows(seqs[k], 11)[0]) for k in refs}
    out, flipped = {}, []
    for k, s in seqs.items():
        own = [fwd[x] for x in refs if x != k][:n_ref]  # never the sequence itself
        ref = np.unique(np.concatenate(own)) if own else np.zeros(0, np.uint64)
        f, r = sc._windows(s, 11)
        if np.intersect1d(r, ref).size > np.intersect1d(f, ref).size:
            s = s.translate(_COMP)[::-1]
            flipped.append(k)
        out[k] = s
    return out, flipped


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def build_tree(seqs: dict[str, str], tools: Tools, work: Path, threads: int = 0,
               max_gap: float = 0.9) -> tuple[str, dict[str, float]]:
    """Align, drop columns with >= ``max_gap`` gaps, infer a tree. Returns (newick, coverage per name)."""
    threads = threads or os.cpu_count() or 1
    names = {f"s{i}": k for i, k in enumerate(seqs)}  # safe tip names
    fa, aln, trim = work / "in.fa", work / "in.aln", work / "trim.aln"
    fa.write_text("".join(f">{t}\n{seqs[k]}\n" for t, k in names.items()))
    if tools.aligner_kind == "famsa":
        subprocess.run([tools.aligner, "-t", str(threads), str(fa), str(aln)], check=True, capture_output=True)
    else:
        with open(aln, "w") as fh:
            subprocess.run([tools.aligner, "--auto", "--thread", str(threads), "--quiet", str(fa)], stdout=fh,
                           check=True)
    a = _read_fasta(aln)
    m = np.array([list(a[t]) for t in names])
    keep = (m == "-").mean(axis=0) < max_gap
    m = m[:, keep]
    trim.write_text("".join(f">{t}\n{''.join(row)}\n" for t, row in zip(names, m)))
    # coverage of the core region (columns where most sequences have data): a fragment of the common
    # region scores ~1 even beside complete mitogenome genes; a record of another region scores ~0
    data = m != "-"
    core = data.mean(axis=0) >= 0.5
    if not core.any():
        core = np.ones(m.shape[1], bool)
    cover = {names[t]: float(row[core].mean()) for t, row in zip(names, data)}
    if tools.tree_kind == "fasttree":
        args = [tools.tree, "-nt", "-gtr", "-gamma", "-quiet"] + (["-fastest"] if len(seqs) > 1000 else [])
        nwk = subprocess.run(args + [str(trim)], check=True, capture_output=True, text=True).stdout
    else:
        pre = work / "iq"
        subprocess.run([tools.tree, "-s", str(trim), "-B", "1000", "-T", str(threads), "--prefix", str(pre), "-quiet",
                        "-redo"], check=True, capture_output=True)
        nwk = (work / "iq.treefile").read_text()
    from io import StringIO

    from Bio import Phylo

    tree = Phylo.read(StringIO(nwk), "newick")
    for c in tree.get_terminals():
        c.name = names[c.name]
    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip(), cover


# ------------------------------------------------------------------ splits and decisions


def splits(newick: str) -> tuple[list[str], list[tuple[int, float]]]:
    """Tip names and the tree's splits as (bitmask of tips on one side, support 0-1)."""
    from io import StringIO

    from Bio import Phylo

    tree = Phylo.read(StringIO(newick), "newick")
    tips = [c.name for c in tree.get_terminals()]
    idx = {n: i for i, n in enumerate(tips)}
    out = []
    for cl in tree.get_nonterminals():
        if cl is tree.root:
            continue
        mask = 0
        for t in cl.get_terminals():
            mask |= 1 << idx[t.name]
        sup = cl.confidence if cl.confidence is not None else 0.0
        out.append((mask, sup / 100 if sup > 1 else sup))
    return tips, out


def neighbours(i: int, n: int, parts: list[tuple[int, float]], usable: int, support: float = 0.9,
               min_others: int = 2) -> tuple[int, float]:
    """The smallest side of a split with support >= ``support`` that holds tip ``i`` and >= ``min_others``
    usable other tips: (mask of those others, support); (0, 0) if none. Root-independent."""
    full = (1 << n) - 1
    best, best_sup, best_size = 0, 0.0, n + 1
    for mask, sup in parts:
        if sup < support:
            continue
        side = mask if mask >> i & 1 else full ^ mask
        others = side & usable & ~(1 << i)
        k = bin(others).count("1")
        size = bin(side).count("1")
        if k >= min_others and size < best_size:
            best, best_sup, best_size = others, sup, size
    return best, best_sup


@dataclass
class Tip:
    name: str
    organism: str
    genus: str
    family: str
    study: bool  # True: the studied taxon; False: background


def judge(tips: list[Tip], parts: list[tuple[int, float]], cover: dict[str, float], study_family: str,
          support: float = 0.9, min_cover: float = 0.5, fast: bool = True) -> tuple[list[dict], list[dict]]:
    """(record flags, taxon notes) for the study tips of one gene tree."""
    n = len(tips)
    usable = sum(1 << i for i, t in enumerate(tips) if cover.get(t.name, 0) >= min_cover)
    members = range(n)
    rows, notes = [], []
    raw: list[dict] = []
    for i, t in enumerate(tips):
        if not t.study or not usable >> i & 1:
            continue
        nb, sup = neighbours(i, n, parts, usable, support)
        if not nb:
            continue
        others = [j for j in members if nb >> j & 1]
        fams = {tips[j].family for j in others}
        mates = [j for j in members if j != i and tips[j].study and t.genus and tips[j].genus == t.genus
                 and usable >> j & 1]
        mates_in = sum(1 for j in mates if nb >> j & 1)
        level = ""
        if study_family not in fams and len(fams) <= 2:
            level, where = "family", "/".join(sorted(fams))
        elif fast and fams == {study_family} and t.genus and mates_in == 0:  # genera: not on 18S / 28S
            ogen = {tips[j].genus for j in others if tips[j].genus}
            if ogen and t.genus not in ogen and all(tips[j].genus for j in others):
                level, where = "genus", "/".join(sorted(ogen))
        if not level:
            continue
        raw.append({"i": i, "level": level, "where": where, "support": round(sup, 2), "mates": len(mates),
                    "mates_in": mates_in, "others": [tips[j].name for j in others],
                    "neighbours": "; ".join(sorted({tips[j].organism for j in others}))[:200]})
    # a genus flagged as a whole (congeners with it, or most of the genus flagged) is a taxonomic note
    by_genus: dict[str, list[dict]] = {}
    for r in raw:
        by_genus.setdefault(tips[r["i"]].genus, []).append(r)
    for g, rs in by_genus.items():
        n_g = sum(1 for j, t in enumerate(tips) if t.study and g and t.genus == g and usable >> j & 1)
        whole = bool(g) and (len(rs) > 0.5 * n_g and n_g > 1)
        for r in rs:
            t = tips[r["i"]]
            if whole or not r["mates"] or r["mates_in"] >= 0.5 * r["mates"]:
                # the genus as a whole, or its only sequence: a question of genus placement / names
                notes.append({"genus": g, "level": r["level"], "where": r["where"], "support": r["support"],
                              "accession": t.name, "organism": t.organism, "others": r["others"],
                              "why": "genus as a whole" if whole or r["mates"] else "only sequence of its genus"})
                continue
            kind = "outside_family" if r["level"] == "family" else "inside_other_genus"
            rows.append({"accession": t.name, "organism": t.organism, "flag": kind, "others": r["others"],
                         "note": f"sits in {r['where']} (support {r['support']}); "
                                 f"{r['mates']} other {t.genus} sequences elsewhere",
                         "support": r["support"], "neighbours": r["neighbours"]})
    return rows, notes


# ------------------------------------------------------------------ data


def _table(g2t_dir: str | Path) -> pd.DataFrame:
    return pd.read_csv(Path(g2t_dir) / "labeled_genes" / "assigned_genes_types_all.csv", dtype=str,
                       keep_default_na=False, low_memory=False)


def study_family(g2t_dir: str | Path) -> str:
    """The commonest Family of the records of a g2t run."""
    a = _table(g2t_dir)
    fam = a["Family"][a["Family"] != ""] if "Family" in a else pd.Series([], dtype=str)
    return str(fam.mode().iloc[0]) if len(fam) else "study"


def marker_sets(inputs: Iterable[str], g2t_dir: str | Path, study: bool,
                family: str = "") -> dict[str, dict[str, tuple[str, Tip]]]:
    """{gene: {accession.version: (sequence, Tip)}} from a g2t run (single-gene records + mitogenome genes)."""
    from g2t.verify import read_sequences

    inputs = list(inputs)
    a = _table(g2t_dir)
    seqs = {acc.split(".")[0]: (acc, s) for acc, _, s in read_sequences(inputs, set(a["ACCESSION"]))}
    out: dict[str, dict[str, tuple[str, Tip]]] = {}

    def tip(row: pd.Series, acc: str) -> Tip:
        fam = row.get("Family", "") or family
        return Tip(acc, row["Organism"], sc.genus_of(row["Organism"], row.get("Genus", "")),
                   family if study else fam, study)

    rows = {r["ACCESSION"]: r for _, r in a.iterrows()}
    for acc, gene in zip(a["ACCESSION"], a["gene_type"]):
        if gene in MARKERS and acc in seqs:
            out.setdefault(gene, {})[seqs[acc][0]] = (seqs[acc][1], tip(rows[acc], seqs[acc][0]))
    mt = set(a.loc[a["gene_type"] == "mtgenome", "ACCESSION"])
    for acc, genes in (sc.mitogenome_genes(inputs, mt) if mt else {}).items():
        for gene, s in genes.items():
            if gene in MARKERS:
                out.setdefault(gene, {})[acc] = (s, tip(rows[acc.split(".")[0]], acc))
    return out


def sample(entries: dict[str, tuple[str, Tip]], per_family: int) -> dict[str, tuple[str, Tip]]:
    """Background: the longest sequence of each species, at most ``per_family`` per family."""
    best: dict[str, tuple[str, str, Tip]] = {}
    for acc, (s, t) in entries.items():
        key = t.organism
        if key not in best or len(s) > len(best[key][1]):
            best[key] = (acc, s, t)
    out: dict[str, tuple[str, Tip]] = {}
    count: dict[str, int] = {}
    for acc, s, t in sorted(best.values(), key=lambda x: -len(x[1])):
        if count.get(t.family, 0) < per_family:
            count[t.family] = count.get(t.family, 0) + 1
            out[acc] = (s, t)
    return out


def representatives(entries: dict[str, tuple[str, Tip]], cap: int, keep: set[str]) -> dict[str, tuple[str, Tip]]:
    """All study sequences if <= ``cap``; else the longest per species plus ``keep`` (distance-flagged)."""
    if len(entries) <= cap:
        return entries
    best: dict[str, str] = {}
    for acc, (s, t) in entries.items():
        if t.organism not in best or len(s) > len(entries[best[t.organism]][0]):
            best[t.organism] = acc
    chosen = set(best.values()) | (keep & set(entries))
    return {a: entries[a] for a in chosen}


def intermixed_genera(rows: list[dict], notes: list[dict]) -> tuple[list[dict], list[dict]]:
    """Genus A inside genus B and B inside A (any marker): the two names are used for one lineage
    (synonyms or unsettled limits, e.g. Ophiodromus / Oxydromus), a taxon note rather than record errors."""
    def gen(r: dict) -> str:
        return sc.genus_of(r["organism"])

    pairs = {(gen(r), w) for r in rows if r["flag"] == "inside_other_genus" for w in r["note"].split(" (")[0]
             .removeprefix("sits in ").split("/")}
    keep = []
    for r in rows:
        where = r["note"].split(" (")[0].removeprefix("sits in ").split("/")
        if r["flag"] == "inside_other_genus" and any((w, gen(r)) in pairs for w in where):
            notes.append({"gene": r["gene"], "genus": gen(r), "organism": r["organism"], "accession": r["accession"],
                          "level": "genus", "where": "/".join(where), "support": r["support"],
                          "why": "genera intermixed: synonyms or unsettled genus limits",
                          "group_distance": r.get("group_distance")})
        else:
            keep.append(r)
    return keep, notes


def related_families(taxon: str, n: int = 6, gene: str = "COI,16S,18S,28S", min_records: int = 20) -> list[str]:
    """Background families for ``taxon``: walk up its NCBI lineage to the nearest ancestor holding at
    least ``n`` other families, and take the ``n`` with most marker records (e.g. Opheliidae ->
    Scolecida: Capitellidae, Maldanidae, Orbiniidae ...)."""
    from Bio import Entrez

    from g2t import download as dl

    dl.setup_entrez()
    tid, name, _, _ = dl.resolve_taxon(taxon)
    info = dl._call(Entrez.efetch, True, db="taxonomy", id=tid, retmode="xml")[0]
    own = {name} | {x["ScientificName"] for x in info.get("LineageEx", []) if x.get("Rank") == "family"}
    for anc in reversed(info.get("LineageEx", [])):
        r = dl._call(Entrez.esearch, True, db="taxonomy", term=f"txid{anc['TaxId']}[Subtree] AND family[Rank]",
                     retmax=500)
        ids = r["IdList"]
        if len(ids) - 1 < n:
            continue
        fams = [x["ScientificName"] for x in dl._call(Entrez.efetch, True, db="taxonomy", id=",".join(ids),
                                                         retmode="xml")]
        counts = []
        for f, fid in zip(fams, ids):
            if f in own:
                continue
            term, _ = dl.build_query(fid, dl.DownloadOptions(gene=gene))
            counts.append((dl.count(term), f))
        chosen = [f for c, f in sorted(counts, reverse=True) if c >= min_records][:n]
        logger.info("background for %s: %s (from %s %s)", name, ", ".join(chosen), anc.get("Rank", ""),
                    anc["ScientificName"])
        return chosen
    return []


def fetch_background(families: list[str], out_dir: str | Path, gene: str = "COI,16S,18S,28S"
                     ) -> list[tuple[list[str], str]]:
    """Download the families' marker records (incremental) and type their genes (g2t steps 1-2)."""
    from g2t import download as dl
    from g2t.classify import classify
    from g2t.extract import extract

    out = Path(out_dir)
    dirs = []
    for f in families:
        res = dl.download(f, str(out / f), dl.DownloadOptions(gene=gene), progress=False)
        dirs.append(res.out_dir)
    g2t_dir = out / "g2t"
    extract(dirs, str(g2t_dir / "gb_metadata"), stream=True, final_name="final.csv")
    classify(str(g2t_dir / "gb_metadata" / "final.csv"), str(g2t_dir / "labeled_genes"))
    return [(dirs, str(g2t_dir))]


def distance_check(rows: list[dict], study: dict[str, tuple[str, Tip]],
                   entries: dict[str, tuple[str, Tip]]) -> list[dict]:
    """Keep a tree flag only if the sequence is also closer (shared 11-mers) to the group it sits in than
    to its nearest congener. A well-supported placement can still be an artefact: Travisia forbesii
    HQ025027 COI sat in Cirratulidae (support 1.0) but is 3.9% from another T. forbesii and >= 31% from
    any cirratulid; Ophelina acuminata HQ025072 28S sat in Ophelia (0.98) but is identical to another
    O. acuminata in a region too conserved to tell the genera apart."""
    kept = []
    for r in rows:
        acc = r["accession"]
        q = sc.kmers(entries[acc][0], 11)
        genus = entries[acc][1].genus
        group_d = min((sc.pdist(q, sc.kmers(entries[o][0], 11), 11) for o in r["others"]), default=1.0)
        own = [sc.pdist(q, sc.kmers(s, 11), 11) for a, (s, t) in study.items()
               if a != acc and genus and t.genus == genus]
        own_d = min(own, default=None)
        r["group_distance"] = round(group_d, 4)
        r["own_genus_distance"] = None if own_d is None else round(own_d, 4)
        if own_d is None or group_d < own_d:
            kept.append(r)
        else:
            logger.info("%s: tree places it in %s, but it is closer to its own genus (%.1f%% vs %.1f%%): "
                        "not flagged", acc, r["note"].split(" (")[0], 100 * own_d, 100 * group_d)
    return kept


def specimen_check(flags: pd.DataFrame, judged: dict[str, str], matrix: Path) -> pd.DataFrame:
    """Add the other markers of the same specimen: misplaced too (specimen misidentified) or placed
    correctly (this sequence is the problem)."""
    if not len(flags) or not matrix.exists():
        flags["other_markers"] = ""
        flags["evidence"] = ""
        return flags
    m = pd.read_csv(matrix, dtype=str, keep_default_na=False, low_memory=False)
    genes = [c for c in m.columns if c in MARKERS or c == "mtgenome"]
    where: dict[str, int] = {}
    for i, r in m.iterrows():
        for g in genes:
            for a in str(r[g]).replace(";", ",").split(","):
                if a.strip():
                    where[a.strip()] = int(str(i))
    other, evidence = [], []
    for acc in flags["accession"]:
        i = where.get(acc)
        if i is None:
            other.append("")
            evidence.append("specimen not in matrix")
            continue
        r = m.iloc[i]
        parts, states = [], []
        for g in genes:
            for a in str(r[g]).replace(";", ",").split(","):
                a = a.strip()
                if a and a != acc:
                    st = judged.get(a, "not judged")
                    parts.append(f"{g} {a} {st}")
                    states.append(st)
        other.append("; ".join(parts))
        if "flagged" in states:
            evidence.append("other markers of the specimen misplaced too: specimen misidentified?")
        elif "ok" in states:
            evidence.append("other markers of the specimen sit correctly: this sequence is wrong "
                            "(contamination or mix-up)")
        else:
            evidence.append("no other marker judged")
    flags["other_markers"] = other
    flags["evidence"] = evidence
    return flags


def phylo_screen(inputs: list[str], g2t_dir: str | Path, background: list[tuple[list[str], str]],
                 family: str = "", markers: Iterable[str] = MARKERS, out_dir: str | Path | None = None,
                 iqtree: bool = False, threads: int = 0, support: float = 0.9, min_cover: float = 0.5,
                 cap: int = 800, bg_per_family: int = 12) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Record flags and taxon notes for every marker with background sequences."""
    tools = Tools.find(iqtree)
    family = family or study_family(g2t_dir)
    study = marker_sets(inputs, g2t_dir, True, family)
    bg: dict[str, dict[str, tuple[str, Tip]]] = {}
    for gb, d in background:
        for gene, e in marker_sets(gb, d, False).items():
            bg.setdefault(gene, {}).update({k: v for k, v in e.items() if v[1].family != family})
    flagged_before: set[str] = set()
    prev = Path(g2t_dir) / "organized_genes" / "screen.csv"
    if prev.exists():
        flagged_before = set(pd.read_csv(prev)["accession"].astype(str))
    out = Path(out_dir) if out_dir else Path(g2t_dir) / "organized_genes" / "phylo"
    out.mkdir(parents=True, exist_ok=True)
    rows, notes, judged = [], [], {}
    for gene in markers:
        if gene not in study or gene not in bg or len(study[gene]) < 3:
            continue
        ing = representatives(study[gene], cap, flagged_before)
        bgs = sample(bg[gene], bg_per_family)
        entries = {**ing, **bgs}
        seqs, _ = orient({k: v[0] for k, v in entries.items()})
        with tempfile.TemporaryDirectory(dir=out) as w:
            nwk, cover = build_tree(seqs, tools, Path(w), threads)
        (out / f"{gene}.tre").write_text(nwk)
        names, parts = splits(nwk)
        tips = [entries[nm][1] for nm in names]
        r, nt = judge(tips, parts, cover, family, support, min_cover, fast=gene in sc.FAST)
        r = distance_check(r, study[gene], entries)
        for x in nt:  # a lone sequence almost identical to another family is a contaminant, not a taxon question
            q = sc.kmers(entries[x["accession"]][0], 11)
            x["group_distance"] = round(min((sc.pdist(q, sc.kmers(entries[o][0], 11), 11) for o in x["others"]),
                                            default=1.0), 4)
            lone = x["why"] == "only sequence of its genus"
            if x["level"] == "family" and lone and x["group_distance"] <= NEAR_FAMILY:
                r.append({"accession": x["accession"], "organism": x["organism"], "flag": "outside_family",
                          "others": x["others"], "support": x["support"], "group_distance": x["group_distance"],
                          "own_genus_distance": None, "neighbours": "",
                          "note": f"sits in {x['where']} (support {x['support']}), {x['group_distance']:.1%} "
                                  "from it; only sequence of its genus"})
        nt = [x for x in nt if x["accession"] not in {y["accession"] for y in r}]
        flagged = {x["accession"] for x in r}
        for t in tips:
            if t.study:
                judged[t.name] = ("flagged" if t.name in flagged else "ok" if cover.get(t.name, 0) >= min_cover
                                  else "not judged")
        for x in r:
            x["gene"] = gene
        for x in nt:
            x["gene"] = gene
        rows += r
        notes += nt
        logger.info("%s: %d %s + %d background sequences; %d flagged, %d taxon notes", gene, len(ing), family,
                    len(bgs), len(r), len(nt))
    rows, notes = intermixed_genera(rows, notes)
    cols = ["gene", "accession", "organism", "flag", "note", "support", "group_distance", "own_genus_distance",
            "neighbours"]
    flags = specimen_check(pd.DataFrame(rows, columns=cols), judged,
                           Path(g2t_dir) / "organized_genes" / "organized_species_voucher.csv")
    flags.to_csv(out / "phylo_flags.csv", index=False)
    tn = pd.DataFrame(notes, columns=["gene", "genus", "organism", "accession", "level", "where", "why", "support",
                                      "group_distance"])
    tn.to_csv(out / "phylo_taxon_notes.csv", index=False)
    return flags, tn
