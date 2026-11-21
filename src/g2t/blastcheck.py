"""g2t-blastcheck: where does BLAST say a sequence belongs, and does that disagree with its name?

For each query sequence the BLAST hits (outfmt 6 with taxids, see ``OUTFMT``) are reduced to a
taxon-aware consensus (a weighted LCA, as in MEGAN / BASTA) and compared with the lineage the record
declares. The result is evidence, not a verdict: g2t-screen's other layers (k-mer, tree) and the
final judgement add to it.

Rules that keep it honest:
  * the query's own record (same accession, version ignored) is dropped from its hits;
  * hits are merged into *independent units*: one species submitted in one batch of consecutive
    accessions counts once, so 200 sequences of a well-sampled genus do not outvote 1 of a rare one,
    and a batch of contaminants cannot vouch for itself;
  * only hits within 5% of the best bitscore and covering >= 50% of the query vote; a rank is
    accepted when >= 80% of the units agree;
  * the identity of the best hit limits how fine a rank may be judged (a 78% COI hit cannot say
    "genus"); below the lowest threshold the result is ``no_close_reference`` (new lineage or sparse
    database), never a conflict;
  * a conflict needs the foreign taxon to beat the declared genus / family by ``MARGIN`` points and to
    rest on >= 2 independent units (or 1 unit from the query's own batch: a batch mix-up).
"""

from __future__ import annotations

import logging
import re
import sqlite3
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

logger = logging.getLogger("g2t.blastcheck")

OUTFMT = "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore staxids sscinames sstart send"
COLS = OUTFMT.split()[1:]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
TOP = 0.95  # votes: bitscore >= TOP * best
MIN_QCOV = 50.0
AGREE = 0.80  # share of units that must agree for a rank
MARGIN = 2.0  # foreign best identity must exceed the declared taxon's by this many points
DISTANT = 85.0  # phylum / class conflicts need a best hit at least this close
SINGLE_PIDENT = 97.0  # one foreign unit is enough if it is this close and far better than the declared taxon
SINGLE_MARGIN = 10.0
BATCH_SPAN = 500  # accessions of one prefix this close are one submission batch

# best-hit identity needed to judge a rank, per marker class
FAST_LIMITS = {"species": 97.0, "genus": 90.0, "family": 82.0, "order": 75.0, "class": 75.0, "phylum": 0.0}
CONSERVED_LIMITS = {"species": 101.0, "genus": 101.0, "family": 99.0, "order": 97.0, "class": 95.0, "phylum": 90.0}
LOWEST = {"fast": 75.0, "conserved": 90.0}  # below this: no_close_reference
CONSERVED = {"18s", "28s", "18-28s", "h3", "h4", "h2a"}


# ------------------------------------------------------------------ taxonomy (ete3's taxa.sqlite, read directly)


class Taxonomy:
    """Lineage by taxid from the NCBI taxonomy database that ete3 keeps in ~/.etetoolkit/taxa.sqlite."""

    def __init__(self, path: str | Path | None = None):
        p = Path(path) if path else Path.home() / ".etetoolkit" / "taxa.sqlite"
        self.ok = p.exists()
        self.db = sqlite3.connect(f"file:{p}?mode=ro", uri=True) if self.ok else None
        self._node: dict[int, tuple[int, str, str]] = {}
        self._lin: dict[int, dict[str, str]] = {}

    def node(self, taxid: int) -> tuple[int, str, str] | None:
        if taxid not in self._node and self.db is not None:
            row = self.db.execute("SELECT parent, rank, spname FROM species WHERE taxid=?", (taxid,)).fetchone()
            self._node[taxid] = (int(row[0] or 0), row[1], row[2]) if row else None  # type: ignore[assignment]
        return self._node.get(taxid)

    def lineage(self, taxid: int) -> dict[str, str]:
        """{rank: name} for the ranks in ``RANKS`` (those the taxon has)."""
        if taxid in self._lin:
            return self._lin[taxid]
        out: dict[str, str] = {}
        t, seen = taxid, set()
        while t and t not in seen:
            seen.add(t)
            n = self.node(t)
            if not n:
                break
            parent, rank, name = n
            if rank in RANKS:
                out[rank] = name
            if parent == t:
                break
            t = parent
        self._lin[taxid] = out
        return out

    def by_name(self, name: str) -> dict[str, str]:
        """Lineage of a taxon given by scientific name ({} if unknown or the name is ambiguous)."""
        if not self.ok or not name:
            return {}
        rows = self.db.execute("SELECT taxid FROM species WHERE spname=?", (name,)).fetchall()  # type: ignore[union-attr]
        if len(rows) != 1:
            return {}
        return self.lineage(int(rows[0][0]))


# ------------------------------------------------------------------ hits


def read_hits(path: str | Path) -> pd.DataFrame:
    with open(path) as fh:
        ncol = len(fh.readline().rstrip("\n").split("\t")) if Path(path).stat().st_size else len(COLS)
    names = COLS[:ncol]  # COLS starts with qseqid; older files have no sstart / send
    df = pd.read_csv(path, sep="\t", header=None, names=names, dtype={"staxids": str, "sseqid": str},
                     keep_default_na=False)
    if "sstart" not in df:
        df["sstart"], df["send"] = 0, 1
    return df


def acc_key(acc: str) -> str:
    """Accession without version, and without a ``|gb|`` style wrapper."""
    a = acc.split("|")[-2] if acc.count("|") >= 2 else acc
    return a.split(".")[0]


def batch_of(acc: str) -> tuple[str, int]:
    """(prefix, number // BATCH_SPAN): consecutive accessions of one prefix are one batch (HM7467xx)."""
    m = re.match(r"([A-Za-z_]+)(\d+)", acc_key(acc))
    return (m.group(1), int(m.group(2)) // BATCH_SPAN) if m else (acc_key(acc), 0)


def same_batch(a: str, b: str) -> bool:
    pa, na = batch_of(a)
    pb, nb = batch_of(b)
    if pa != pb:
        return False
    ma, mb = re.search(r"\d+", acc_key(a)), re.search(r"\d+", acc_key(b))
    return bool(ma and mb and abs(int(ma.group()) - int(mb.group())) <= BATCH_SPAN)


@dataclass
class Result:
    accession: str
    lca_rank: str = ""  # finest rank at which the units agree
    lca_taxon: str = ""
    best_pident: float = 0.0
    n_units: int = 0
    gate: str = "no_hit"  # actionable | no_close_reference | no_hit
    conflict_rank: str = ""  # highest rank where the consensus differs from the declared lineage
    own_best: float = 0.0  # best identity to a hit of the declared genus
    own_species_best: float = 0.0  # best identity to a hit of the declared species
    other_best: float = 0.0
    same_batch_hit: bool = False
    agrees: bool = False  # consensus matches the declared lineage (positive evidence)
    confirmed_rank: str = ""  # finest rank at which the consensus equals the declared lineage
    declared_rank: str = ""  # finest of species / genus / family the record declares
    resolved: bool = False  # confirmed_rank reaches the declared rank: the name itself is confirmed
    reverse: bool = False  # best hit is on the minus strand: the record is probably reverse-complemented
    note: str = ""
    units: list = field(default_factory=list, repr=False)


OPEN_EPITHETS = {"sp", "cf", "aff", "nr", "gen", "environmental", "uncultured"}


def species_of(organism: str) -> str:
    """'Genus species' of an organism name; '' for 'Genus sp.', 'cf.', 'aff.', 'environmental sample' etc."""
    w = str(organism).split()
    if len(w) < 2 or not w[0][:1].isupper() or w[1] in OPEN_EPITHETS or not re.fullmatch(r"[a-z][a-z-]+", w[1]):
        return ""
    return f"{w[0]} {w[1]}"


def limits_for(gene: str) -> tuple[dict[str, float], float]:
    cons = gene.lower() in CONSERVED
    return (CONSERVED_LIMITS, LOWEST["conserved"]) if cons else (FAST_LIMITS, LOWEST["fast"])


def check_query(q: str, hits: pd.DataFrame, declared: dict[str, str], gene: str, tax: Taxonomy) -> Result:
    """Consensus of the hits of one query against its declared lineage ({rank: name})."""
    res = Result(accession=q)
    h = hits[(hits["sseqid"].map(acc_key) != acc_key(q)) & (hits["qcovs"] >= MIN_QCOV)]
    if h.empty:
        return res
    limits, lowest = limits_for(gene)
    best = h["bitscore"].max()
    # only the hits that vote (bitscore within TOP of the best) count for identity: a shorter, lower-scoring
    # 99% alignment must not unlock a fine rank that the voting references do not support
    voters = h[h["bitscore"] >= TOP * best]
    res.best_pident = float(voters["pident"].max())
    top = h.loc[h["bitscore"].idxmax()]
    res.reverse = bool("sstart" in top and top["sstart"] > top["send"])
    if res.best_pident < lowest:
        res.gate = "no_close_reference"
        return res
    # independent units: (species taxid, batch); the unit keeps its best identity
    units: dict[tuple, dict] = {}
    own_genus = declared.get("genus", "")
    own_family = declared.get("family", "")
    own_species = declared.get("species", "")
    for r in h.itertuples():
        tid = str(r.staxids).split(";")[0]
        if not tid.isdigit():
            continue
        lin = tax.lineage(int(tid))
        if not lin:
            continue
        sb = same_batch(q, r.sseqid)
        if own_genus and lin.get("genus") == own_genus:
            res.own_best = max(res.own_best, float(r.pident))
        if own_family and lin.get("family") == own_family:
            res.own_best = max(res.own_best, float(r.pident)) if not own_genus else res.own_best
        if own_species and lin.get("species") == own_species:
            res.own_species_best = max(res.own_species_best, float(r.pident))
        if r.bitscore < TOP * best:
            continue
        key = (tid, batch_of(r.sseqid))
        u = units.setdefault(key, {"lin": lin, "pident": 0.0, "same_batch": False, "acc": r.sseqid})
        u["pident"] = max(u["pident"], float(r.pident))
        u["same_batch"] = u["same_batch"] or sb
    res.n_units = len(units)
    if not units:
        res.note = "taxids_unknown: hit taxids missing from the local taxonomy (update ete3: ete3 ncbiupdate)"
        return res
    res.units = list(units.values())
    res.gate = "actionable"
    # finest accepted rank, and the finest rank the identity allows
    allowed = [rk for rk in RANKS if res.best_pident >= limits[rk]]
    finest_ok = RANKS.index(allowed[-1]) if allowed else 0
    consensus: dict[str, str] = {}
    for rk in RANKS[: finest_ok + 1]:
        votes: dict[str, int] = defaultdict(int)
        for u in units.values():
            name = u["lin"].get(rk)
            if rk == "species" and name and not species_of(name):
                name = ""  # 'Nereis aff. zonata', 'Nereididae sp. X': an open name is no species
            if name:  # a hit without this rank (no 'order' in annelids, 'Polychaeta sp.') abstains
                votes[u["lin"][rk]] += 1
        if not votes:
            continue
        top, n = max(votes.items(), key=lambda kv: kv[1])
        if n / sum(votes.values()) >= AGREE:
            consensus[rk] = top
            res.lca_rank, res.lca_taxon = rk, top
        elif rk in ("genus", "species"):
            break  # fine ranks disagree: stop here, coarser ranks stay
        else:
            break
    # conflict with the declared lineage: highest rank where they differ
    for rk in RANKS:
        if rk in consensus and declared.get(rk) and consensus[rk] != declared[rk]:
            res.conflict_rank = rk
            break
    if res.conflict_rank:
        foreign = [u for u in units.values() if u["lin"].get(res.conflict_rank) == consensus[res.conflict_rank]]
        res.other_best = max(u["pident"] for u in foreign)
        res.same_batch_hit = any(u["same_batch"] for u in foreign)
        # species-level conflict (same genus, other species): compare with the declared *species* only,
        # not the whole genus, which contains the foreign species itself
        own_cmp = res.own_species_best if res.conflict_rank == "species" else res.own_best
        ok_margin = res.other_best - own_cmp >= MARGIN
        ok_units = len(foreign) >= 2 or res.same_batch_hit
        single = (not ok_units and res.other_best >= SINGLE_PIDENT and res.other_best - own_cmp >= SINGLE_MARGIN)
        ok_units = ok_units or single
        if res.conflict_rank in ("phylum", "class") and res.best_pident < DISTANT:
            res.note = f"distant_hit_only: best hit below {DISTANT:.0f}%, a new lineage cannot be told from a stranger"
            res.conflict_rank = ""
        elif not (ok_margin and ok_units):
            res.note = "discordant but not independent / not clearly better than declared taxon"
            res.conflict_rank = ""
        elif res.own_best == 0 and res.conflict_rank in ("genus", "family"):
            # the declared taxon is absent from the top hits: still a conflict, but a weaker one
            res.note = "declared taxon absent from top hits"
        elif res.conflict_rank == "species":
            res.note = f"species_conflict: best hit {consensus['species']} at {res.other_best:.1f}%"
        elif single and len(foreign) < 2 and not res.same_batch_hit:
            res.note = "single_reference: one foreign hit, very close"
    elif consensus:
        # positive evidence: consensus reaches family or finer and agrees with the declared lineage
        both = [rk for rk in consensus if declared.get(rk)]
        res.agrees = any(RANKS.index(rk) >= RANKS.index("family") for rk in both) and not res.note
        agreeing = [rk for rk in both if consensus[rk] == declared[rk]]
        res.confirmed_rank = agreeing[-1] if agreeing else ""
        # a rank this marker can never resolve (conserved rRNA / histone at species level) is not demanded
        fine = [rk for rk in ("family", "genus", "species") if declared.get(rk) and limits[rk] <= 100]
        res.declared_rank = fine[-1] if fine else ""
        res.resolved = bool(res.agrees and res.declared_rank and res.confirmed_rank == res.declared_rank)
    return res


def check_all(hits: pd.DataFrame, declared: dict[str, dict[str, str]], genes: dict[str, str],
              tax: Taxonomy | None = None) -> list[Result]:
    """``declared``: accession -> {rank: name}; ``genes``: accession -> gene type."""
    tax = tax or Taxonomy()
    out = []
    for q, g in hits.groupby("qseqid", sort=False):
        out.append(check_query(q, g, declared.get(q, declared.get(acc_key(q), {})), genes.get(q, "coi"), tax))
    return out


# ------------------------------------------------------------------ running BLAST


def run_blast(fasta: str | Path, db: str | Path, out: str | Path, threads: int = 16, blastn: str = "blastn",
              max_targets: int = 30) -> Path:
    """megablast of ``fasta`` against a nucleotide database; writes ``out`` (outfmt ``OUTFMT``)."""
    db = Path(db)
    env = {"BLASTDB": str(db.parent)}
    cmd = [blastn, "-task", "megablast", "-db", db.name, "-query", str(fasta), "-num_threads", str(threads),
           "-max_target_seqs", str(max_targets), "-outfmt", OUTFMT, "-out", str(out)]
    logger.info("BLAST: %s", " ".join(cmd))
    import os

    subprocess.run(cmd, check=True, env={**os.environ, **env})
    return Path(out)


# ------------------------------------------------------------------ whole g2t run


def run_blastcheck(inputs, g2t_dir, db, threads: int = 16, blastn: str = "blastn", genes=None,
                   tax: Taxonomy | None = None) -> pd.DataFrame:
    """BLAST every marker of a g2t run against ``db`` and compare with the declared Family / Genus."""
    from g2t import phylo

    d = Path(g2t_dir)
    work = d / "organized_genes" / "blastcheck"
    work.mkdir(parents=True, exist_ok=True)
    tax = tax or Taxonomy()
    a = phylo._table(d)
    rows = {r["ACCESSION"]: r for _, r in a.iterrows()}
    sets = phylo.marker_sets(list(inputs), d, True, "")
    out: list[dict] = []
    for g, entries in sets.items():
        if genes and g not in genes or not entries:
            continue
        fa, hf = work / f"{g}.fa", work / f"{g}.hits.tsv"
        with open(fa, "w") as fh:
            for acc, (seq, _tip) in entries.items():
                fh.write(f">{acc}\n{seq}\n")
        run_blast(fa, db, hf, threads=threads, blastn=blastn)
        if hf.stat().st_size == 0:
            continue
        hits = read_hits(hf)
        declared, org = {}, {}
        for q in hits["qseqid"].unique():
            r = rows.get(acc_key(q))
            if r is None:
                continue
            org[q] = r["Organism"]
            lin = tax.by_name(r["Family"]) or tax.by_name(r["Genus"])
            declared[q] = {**lin, **({"family": r["Family"]} if r["Family"] else {}),
                           **({"genus": r["Genus"]} if r["Genus"] else {})}
            sp = species_of(r["Organism"])
            if sp:
                declared[q]["species"] = sp
        for x in check_all(hits, declared, {q: g for q in declared}, tax):
            out.append({"accession": x.accession, "gene": g, "organism": org.get(x.accession, ""), "gate": x.gate,
                        "lca_rank": x.lca_rank, "lca_taxon": x.lca_taxon, "best_pident": x.best_pident,
                        "n_units": x.n_units, "conflict_rank": x.conflict_rank, "own_best": x.own_best,
                        "other_best": x.other_best, "same_batch_hit": x.same_batch_hit, "agrees": x.agrees,
                        "reverse": x.reverse, "note": x.note, "confirmed_rank": x.confirmed_rank,
                        "declared_rank": x.declared_rank, "resolved": x.resolved})
    return pd.DataFrame(out)
