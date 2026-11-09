from __future__ import annotations

import argparse
import contextlib
import importlib.util
import logging
import os
import re
import sys
import warnings
from dataclasses import dataclass, field
from functools import cache, lru_cache
from typing import Any

import pandas as pd

from g2t.genes import MITO_TYPES, NUCLEAR_TYPES, PLASTID_TYPES
from g2t.utils import StepResult, UnsupportedFormatError, length_in_range, parse_interval
from g2t.utils import read_table as utils_read_table


# =========================
# Logging configuration
# =========================
class CustomFormatter(logging.Formatter):
    def format(self, record):
        if record.levelno == logging.DEBUG:
            record.levelname = "CHECKING"
        return super().format(record)


def setup_logging(log_file: str = "classify.log") -> logging.Handler:
    """Write the full (DEBUG) log of this run to ``log_file``. The console is left to the caller
    (``g2t-classify`` configures it once; adding a second console handler printed every line twice)."""
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(CustomFormatter("%(asctime)s - %(levelname)s - %(message)s"))
    pkg_logger = logging.getLogger("g2t")
    pkg_logger.setLevel(logging.DEBUG)
    pkg_logger.addHandler(handler)
    return handler


logger = logging.getLogger(__name__)


# =========================
# Gene dictionary: the packaged data/gene_dict.yaml is the single source of synonyms
# =========================
def _yaml_path() -> str:
    try:
        from importlib.resources import files

        return str(files("g2t") / "data" / "gene_dict.yaml")
    except Exception:  # noqa: BLE001 - zipped or unusual installs
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gene_dict.yaml")


def _load_gene_dict_yaml(yaml_path: str | None = None) -> dict[str, list[str]]:
    """Read a gene dictionary in the gene_dict.yaml layout (``gene_types: {gene: {synonyms: [...]}}``).

    Raises on a missing or unreadable file: classifying with an empty or partial dictionary
    would silently leave records unmatched."""
    import yaml

    path = yaml_path or _yaml_path()
    with open(path, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    result: dict[str, list[str]] = {}
    for k, v in (data.get("gene_types") or {}).items():
        synonyms = v.get("synonyms", []) if isinstance(v, dict) else (v if isinstance(v, list) else [])
        result[str(k).lower()] = sorted({s.lower() for s in synonyms if isinstance(s, str)})
    if not result:
        raise ValueError(f"No gene_types in gene dictionary {path}")
    return result


@lru_cache(maxsize=1)
def _packaged_gene_dict() -> tuple[tuple[str, tuple[str, ...]], ...]:
    return tuple((k, tuple(v)) for k, v in _load_gene_dict_yaml().items())


def default_gene_dict() -> dict[str, list[str]]:
    """A fresh, caller-owned copy of the packaged gene dictionary (gene -> lower-case synonyms)."""
    return {k: list(v) for k, v in _packaged_gene_dict()}


# Default gene sets; a run uses the copies in MatchConfig (also the defaults of init_gene_groups)
_DEFAULT_MITO_GENES = set(MITO_TYPES)
_DEFAULT_NUCLEAR_GENES = set(NUCLEAR_TYPES)
_DEFAULT_PLASTID_GENES = set(PLASTID_TYPES)


# =========================
# Keyword patterns (pure memoisation: a phrase always compiles to the same pattern)
# =========================
@cache
def _phrase_pattern(key: str) -> re.Pattern:
    escaped = re.escape(key)
    left = r"\b" if re.match(r"\w", key) else r"(?:(?<=\s)|^)"
    right = r"\b" if (key and re.match(r"\w", key[-1])) else r"(?=\s|$)"
    return re.compile(left + escaped + right, re.IGNORECASE)


def init_gene_groups(
    mtgenes_list: list[str],
    ntgenes_list: list[str],
    cpgenes_list: list[str] | None = None,
    gene_dict: dict[str, list[str]] | None = None,
) -> tuple[set[str], ...]:
    """Normalise the gene lists of a run. Returns (mito_genes, nuclear_genes), plus plastid_genes when
    cpgenes_list is given. Genes without synonyms in ``gene_dict`` (default: the packaged dictionary)
    are logged; they stay in the sets but can never match. Nothing global is changed."""
    known = default_gene_dict() if gene_dict is None else gene_dict
    mito_genes = {g.strip().lower() for g in mtgenes_list}
    nuclear_genes = {g.strip().lower() for g in ntgenes_list}
    plastid_genes = {g.strip().lower() for g in cpgenes_list or []}
    for label, genes in (("Mitochondrial", mito_genes), ("Nuclear", nuclear_genes), ("Plastid", plastid_genes)):
        for gene in sorted(genes):
            if not known.get(gene):
                logger.warning(f"{label} gene '{gene}' has no synonyms in the gene dictionary; it cannot match")

    logger.info(f"Mitochondrial genes: {', '.join(sorted(mito_genes))}")
    logger.info(f"Nuclear genes: {', '.join(sorted(nuclear_genes))}")
    if cpgenes_list is None:
        return mito_genes, nuclear_genes
    logger.info(f"Plastid genes: {', '.join(sorted(plastid_genes))}")
    return mito_genes, nuclear_genes, plastid_genes


# =========================
# Configuration dataclass (holds the gene sets of a run)
# =========================
@dataclass
class MatchConfig:
    if_remove_matched: str = "no"

    length2_mtgenes: str = "150:50000"
    length2_ntgenes: str = "150:50000"
    length2_cpgenes: str = "150:300000"  # plastomes are ~120-220 kb; plastid records bypass the global filter
    global_length_range: str = "150:50000"

    # per-type length range; a type not listed has none ("none")
    gene_length_ranges: dict[str, str] = field(
        default_factory=lambda: {
            "mtgenome": "3000:500000",
            "plastome": "20000:300000",
            "18-28s": "5000:15000",
        }
    )

    add_assignment_reasons: bool = True
    record_conflict_reasons: bool = True

    # gene sets of this run
    mito_genes: set[str] = field(default_factory=lambda: _DEFAULT_MITO_GENES.copy())
    nuclear_genes: set[str] = field(default_factory=lambda: _DEFAULT_NUCLEAR_GENES.copy())
    plastid_genes: set[str] = field(default_factory=lambda: _DEFAULT_PLASTID_GENES.copy())

    # gene -> lower-case synonyms; each config owns its copy (see load_path_dict to add entries)
    gene_dict: dict[str, list[str]] = field(default_factory=default_gene_dict)

    def get_gene_length_range(self, gene: str) -> str:
        return self.gene_length_ranges.get(gene.lower(), "none")


# =========================
# Load external gene dictionary
# =========================
def read_gene_dict_file(path_dict: str) -> dict[str, list[str]]:
    """Read a user gene dictionary (.txt ``gene: syn 1, syn 2`` lines, .csv with ``gene,synonyms``,
    or .py defining ``gene_dict``). Raises on a missing, unreadable or empty file."""
    if not os.path.exists(path_dict):
        raise FileNotFoundError(f"Gene dictionary file not found: {path_dict}")
    ext = os.path.splitext(path_dict)[1].lower()
    new_dict: dict[str, list[str]] = {}
    if ext == ".py":
        spec = importlib.util.spec_from_file_location("ext_dict", path_dict)
        if spec is None or spec.loader is None:
            raise ImportError(f"cannot load {path_dict}")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        if not hasattr(mod, "gene_dict"):
            raise ValueError(f"{path_dict}: Python file has no 'gene_dict' variable")
        new_dict = {str(k).strip().lower(): [str(x).lower() for x in v] for k, v in mod.gene_dict.items()}
    elif ext == ".txt":
        with open(path_dict, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and ":" in line:
                    gene, syns = line.split(":", 1)
                    new_dict[gene.strip().lower()] = [s.strip().lower() for s in syns.split(",") if s.strip()]
    elif ext == ".csv":
        df = pd.read_csv(path_dict, dtype=str)
        if not {"gene", "synonyms"}.issubset(df.columns):
            raise ValueError(f"{path_dict}: CSV must contain 'gene' and 'synonyms' columns")
        for row in df.itertuples(index=False):
            gene = row.gene.strip().lower()
            new_dict[gene] = [s.strip().lower() for s in str(row.synonyms).split(",") if s.strip()]
    else:
        raise ValueError(f"Unsupported gene dictionary format: {ext} (use .txt, .csv or .py)")
    if not new_dict:
        raise ValueError(f"No gene entries in {path_dict}")
    return new_dict


def load_path_dict(path_dict: str, base: dict[str, list[str]] | None = None) -> dict[str, list[str]]:
    """Return ``base`` (default: the packaged dictionary) merged with the entries of ``path_dict``.

    Synonyms of a known gene are added to it; new genes are added. ``base`` is not modified.
    Use the result as ``MatchConfig(gene_dict=...)``."""
    merged = {k: list(v) for k, v in (default_gene_dict() if base is None else base).items()}
    new_dict = read_gene_dict_file(path_dict)
    for gene, synonyms in new_dict.items():
        merged[gene] = sorted(set(merged.get(gene, [])) | set(synonyms))
    logger.info(f"External gene dictionary loaded ({len(new_dict)} genes)")
    return merged


def smart_match(phrase: str, text: str) -> bool:
    """Whole-word, case-insensitive phrase match."""
    return _phrase_pattern(phrase.lower()).search(text) is not None


def any_match(word_list: list[str], text: str) -> bool:
    return any(smart_match(w, text) for w in word_list)


def clean_length(length_str: str) -> int:
    """'550 bp' / '550' / 550.0 / '1,234 bp' -> int; anything else -> 0.
    (Stripping all non-digits turned 550.0 into 5500.)"""
    m = re.search(r"\d[\d,]*(?:\.\d+)?", str(length_str))
    if not m:
        return 0
    try:
        return int(float(m.group(0).replace(",", "")))
    except ValueError:
        return 0


# Words that look like gene keywords but are not marker genes.
NEGATIVE_CONTEXT: dict[str, re.Pattern] = {
    "16s": re.compile(
        r"16s\s+r?rna\s+(?:\((?:[^()]|\([^()]*\))*\)[\s-]*)?"
        r"(?:methyl|methylase|methyltransferase)",
        re.I,
    ),
    "12s": re.compile(r"12s\s+r?rna\s+(?:methyl|methylase|methyltransferase)", re.I),
    "h3": re.compile(
        r"(?:histone\s+)?h3\s*(?:k\d+|lysine|\S*(?:demethylase|methyltransferase|acetyltransferase|"
        r"deacetylase|kinase|chaperone))",
        re.I,
    ),
    "ef-1": re.compile(r"elongation\s+factor[\s-]*1[\s-]*(?:beta|gamma|delta|b|g|d)\b", re.I),
}
_SPECIMEN_TOKEN = re.compile(
    r"\b(?:isolate|voucher|strain|clone|specimen|haplotype|culture|individual|sample)"
    r"\s+(?:voucher\s+)?[^\s,;]+",
    re.I,
)


def clean_definition(definition: Any) -> str:
    """Lower-cased DEFINITION for keyword matching: specimen identifiers ('isolate CO2', 'voucher COI-12')
    and the English word 'its' are removed first, so they are not read as gene names."""
    text = str(definition or "")
    text = _SPECIMEN_TOKEN.sub(" ", text)
    text = re.sub(r"\bits\b", " ", text)  # lower-case possessive only; 'ITS' is kept
    return re.sub(r"\s+", " ", text).lower().strip()


def _negative_hit(gene: str, definition: str) -> bool:
    pat = NEGATIVE_CONTEXT.get(gene)
    return bool(pat and pat.search(definition))


def filter_dataframe(
    df: pd.DataFrame, moleculetype: str, mol_type: str, length_filter: str, organelle_filter: str
) -> pd.DataFrame:
    if moleculetype.lower() != "all":
        df = df[df["MoleculeType"].str.lower() == moleculetype.lower()]
    if mol_type.lower() != "all":
        allowed = [x.strip().lower() for x in mol_type.split(",")]
        df = df[df["mol_type"].str.lower().isin(allowed)]
    if length_filter.lower() != "all":
        lo, hi = parse_interval(length_filter)
        if lo is not None:
            df = df[df["Length"] >= lo]
        if hi is not None:
            df = df[df["Length"] <= hi]
    if organelle_filter.lower() != "all":
        df = df[df["organelle"].str.lower() == organelle_filter.lower()]
    return df


FILTERED_FILE = "filtered_records.csv"
STATUS_FILE = "record_status.csv"
_FILTER_COLS = ["input_row", "ACCESSION", "LocusID", "Organism", "Length", "Definition", "step", "reason"]


def _write_filtered(output_dir: str, rows: list[dict], append: bool) -> None:
    """Write records removed by a filter, with the reason, so nothing disappears silently."""
    path = os.path.join(output_dir, FILTERED_FILE)
    df = pd.DataFrame(rows, columns=_FILTER_COLS)
    if append and os.path.exists(path):
        df = pd.concat([pd.read_csv(path, dtype=str), df.astype(str)], ignore_index=True)
    df.to_csv(path, index=False)


def _filtered_rows(removed: pd.DataFrame, step: str, reason) -> list[dict]:
    out = []
    for _, r in removed.iterrows():
        out.append(
            {
                "input_row": r.get("g2t_input_row", ""),
                "ACCESSION": r.get("ACCESSION", ""),
                "LocusID": r.get("LocusID", ""),
                "Organism": r.get("Organism", r.get("organism", "")),
                "Length": r.get("g2t_length_raw", "") if pd.isna(r.get("Length")) else int(r.get("Length")),
                "Definition": r.get("Definition", ""),
                "step": step,
                "reason": reason(r) if callable(reason) else reason,
            }
        )
    return out


def _filter_dataframe_logged(
    df: pd.DataFrame, moleculetype: str, mol_type: str, length_filter: str, organelle_filter: str, removed: list[dict]
) -> pd.DataFrame:
    """filter_dataframe, but every removed row is appended to `removed` with its reason."""

    def apply(d, keep, reason):
        removed.extend(_filtered_rows(d[~keep], "classify round 1", reason))
        return d[keep]

    if moleculetype.lower() != "all":
        df = apply(
            df,
            df["MoleculeType"].str.lower() == moleculetype.lower(),
            lambda r: f"MoleculeType '{r.get('MoleculeType', '')}' is not '{moleculetype}'",
        )
    if mol_type.lower() != "all":
        allowed = [x.strip().lower() for x in mol_type.split(",")]
        df = apply(
            df,
            df["mol_type"].str.lower().isin(allowed),
            lambda r: f"mol_type '{r.get('mol_type', '')}' not in {mol_type}",
        )
    if length_filter.lower() != "all":
        lo, hi = parse_interval(length_filter)
        keep = pd.Series(True, index=df.index)
        if lo is not None:
            keep &= df["Length"] >= lo
        if hi is not None:
            keep &= df["Length"] <= hi
        df = apply(df, keep, lambda r: f"length {int(r['Length'])} bp outside --length {length_filter}")
    if organelle_filter.lower() != "all":
        df = apply(
            df,
            df["organelle"].str.lower() == organelle_filter.lower(),
            lambda r: f"organelle '{r.get('organelle', '')}' is not '{organelle_filter}'",
        )
    return df


CLASSIFY_OUTPUTS = [
    "assigned_genes_type.csv",
    "assigned_genes_types2.csv",
    "assigned_genes_types_all.csv",
    "unmatched_sequences.csv",
    "unmatched_sequences2.csv",
    "unmatched_categorized.csv",
    "unmatched_report.txt",
    "conflicted_genes.csv",
    "multiple_matched.csv",
    "multiple_matched2.csv",
    FILTERED_FILE,
    STATUS_FILE,
]


def _clean_outputs(output_dir: str, extra: list[str] | None = None) -> None:
    """Remove the products of a previous run, so stale files are never read back."""
    for name in CLASSIFY_OUTPUTS + list(extra or []):
        path = os.path.join(output_dir, name)
        if os.path.isfile(path):
            os.remove(path)


def _write_all(output_dir: str, prematch_output: str) -> str:
    """assigned_genes_types_all.csv = round 1 (+ round 2 if it ran); always written."""
    frames = []
    for name in (prematch_output, "assigned_genes_types2.csv"):
        path = os.path.join(output_dir, name)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            with contextlib.suppress(pd.errors.EmptyDataError):
                frames.append(pd.read_csv(path, dtype=str))
    merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=DIAG_COLUMNS + ["LocusID"])
    if "LocusID" in merged.columns:
        merged = merged.drop_duplicates(subset=["LocusID"])
    merged = reorder_columns(merged)
    all_file = os.path.join(output_dir, "assigned_genes_types_all.csv")
    merged.to_csv(all_file, index=False)
    return all_file


DIAG_COLUMNS = ["gene_type", "Conflict", "match_source", "Original_match", "Conflict_reason", "Assignment_reason"]


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    front = [c for c in DIAG_COLUMNS if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest]


# =========================
# Shared gene matching loop
# =========================
def _match_gene_loop(
    definition: str,
    length_val: int,
    config: MatchConfig,
    gene_set: set[str],
    skip_genes: set[str],
    priority_genes: list[str],
    global_range: str,
    extra_priority_check=None,
) -> tuple[set[str], list[str]]:
    """Shared matching loop for both mitochondrial and nuclear gene sets."""
    matched: set[str] = set()
    notes: list[str] = []

    if not length_in_range(length_val, global_range):
        notes.append(f"Length {length_val} outside range {global_range}")
        return matched, notes

    # Priority genes
    for pg in priority_genes:
        if pg in skip_genes or pg not in gene_set:
            continue
        if any_match(config.gene_dict.get(pg, []), definition) and not _negative_hit(pg, definition):
            range_str = config.get_gene_length_range(pg)
            if length_in_range(length_val, range_str):
                # Extra check for mtgenome
                mt_min = _mtgenome_min(config)
                if pg == "mtgenome" and mt_min is not None and length_val < mt_min:
                    notes.append(f"mtgenome: keyword matched but Length={length_val}<{mt_min}")
                    continue
                matched.add(pg)
                notes.append(f"{pg}: priority match hit")
                return matched, notes
            else:
                notes.append(f"{pg}: keyword matched but length {length_val} outside range {range_str}")

    # Remaining genes
    for gene in sorted(gene_set):
        if gene in skip_genes or gene in priority_genes:
            continue
        keywords = config.gene_dict.get(gene, [])
        if not keywords:
            continue
        if not any_match(keywords, definition):
            continue
        if _negative_hit(gene, definition):
            notes.append(f"{gene}: keyword matched but excluded by context (not the marker gene)")
            continue

        range_str = config.get_gene_length_range(gene)
        if not length_in_range(length_val, range_str):
            notes.append(f"{gene}: keyword matched but length {length_val} outside range {range_str}")
            continue

        matched.add(gene)
        if range_str.lower() != "none":
            notes.append(f"{gene}: keyword matched, length {length_val} in range {range_str}")
        else:
            notes.append(f"{gene}: keyword matched (no length filter)")

        if config.if_remove_matched.lower() == "yes":
            break

    return matched, notes


def match_mito_genes(
    definition: str, length_val: int, molecule_type: str, config: MatchConfig
) -> tuple[set[str], list[str]]:
    """Mitochondrial gene matching: mtgenome priority, then other mito genes."""
    return _match_gene_loop(
        definition,
        length_val,
        config,
        gene_set=config.mito_genes,
        skip_genes=set(),
        priority_genes=["mtgenome"],
        global_range=config.length2_mtgenes,
    )


# rRNA genes flanking the ITS; an ITS record that names one also contains part of it
_ITS_FLANKS = {
    "18s": re.compile(r"\b18s (?:ribosomal rna|rrna)|small subunit ribosomal rna"),
    "28s": re.compile(r"\b2[568]s (?:ribosomal rna|rrna)|large subunit ribosomal rna"),
}


def match_nuclear_genes(definition: str, length_val: int, config: MatchConfig) -> tuple[set[str], list[str]]:
    """Nuclear gene matching: 18-28s priority, then its1-its2, then other nuclear genes.

    An its1-its2 record that also covers part of 18S and/or 28S ("small subunit ribosomal RNA gene and
    internal transcribed spacer 1") is placed in those columns too."""
    matched, notes = _match_gene_loop(
        definition,
        length_val,
        config,
        gene_set=config.nuclear_genes,
        skip_genes=set(),
        priority_genes=["18-28s", "its1-its2"],
        global_range=config.length2_ntgenes,
    )
    if matched == {"its1-its2"}:
        for gene, pat in _ITS_FLANKS.items():
            if gene in config.nuclear_genes and pat.search(definition.lower()):
                matched.add(gene)
                notes.append(f"{gene}: record also covers part of {gene.upper()} (multi-region rDNA)")
    return matched, notes


def match_plastid_genes(definition: str, length_val: int, config: MatchConfig) -> tuple[set[str], list[str]]:
    """Plastid gene matching: plastome first, then the intergenic-spacer markers (psbA-trnH, trnL-trnF,
    rpl32-trnL, ...; their names contain gene names, so they must win over the genes), then the other genes."""
    spacers = ["psba-trnh", "trnl-trnf", "rpl32-trnl", "trns-trng", "atpf-atph", "petb-petd", "psbe-petl", "rbcl-atpb"]
    matched, notes = _match_gene_loop(
        definition,
        length_val,
        config,
        gene_set=config.plastid_genes,
        skip_genes=set(),
        priority_genes=["plastome", *spacers],
        global_range=config.length2_cpgenes,
    )
    # a spacer record may also name a gene ("trnL gene, partial sequence; and maturase K (matK) gene", Angraecum
    # KF672272): the genes are matched too (their synonyms never match a spacer name alone)
    if matched and matched <= set(spacers):
        more, more_notes = _match_gene_loop(
            definition,
            length_val,
            config,
            gene_set=config.plastid_genes,
            skip_genes={"plastome", *spacers},
            priority_genes=[],
            global_range=config.length2_cpgenes,
        )
        matched |= more
        notes += more_notes
    return matched, notes


def _record_group(organelle: str, definition: str) -> str:
    """'mito' when /organelle or the DEFINITION says mitochondrion/mitochondrial; 'plastid' when /organelle is
    plastid/chloroplast or the DEFINITION mentions chloroplast/plastid (not 'nuclear gene for chloroplast
    product'); else 'nuclear'."""
    d = definition.lower()
    if "mitochondrion" in organelle or "mitochondri" in d:
        return "mito"
    if "plastid" in organelle or "chloroplast" in organelle:
        return "plastid"
    if ("chloroplast" in d or "plastid" in d) and "nuclear gene for" not in d:
        return "plastid"
    return "nuclear"


def _mtgenome_min(config: MatchConfig) -> int | None:
    """Configured lower bound of the mtgenome length range; None means no bound."""
    return parse_interval(config.get_gene_length_range("mtgenome"))[0]


def first_round_match(row: dict[str, Any], config: MatchConfig) -> dict[str, Any]:
    definition = clean_definition(row.get("Definition", ""))
    organelle = str(row.get("organelle", "")).lower().strip()
    length_val = clean_length(row.get("Length", "0"))

    group = _record_group(organelle, definition)
    logger.info(f"Round 1 LocusID={row.get('LocusID', 'NA')}, group={group}")
    logger.debug(f"  Definition: {definition[:120]}")

    if group == "mito":
        molecule_type = str(row.get("MoleculeType", "")).lower().strip()
        matched, notes = match_mito_genes(definition, length_val, molecule_type, config)
    elif group == "plastid":
        matched, notes = match_plastid_genes(definition, length_val, config)
    else:
        matched, notes = match_nuclear_genes(definition, length_val, config)

    gene_type = ",".join(sorted(matched))
    match_source = "first_match" if matched else ""
    return _build_result(gene_type, match_source, notes)


# =========================
# Round 1.5: Conflict check
# =========================
def _check_ribosomal_conflict(gene: str, definition: str, length_val: int) -> dict | None:
    if gene == "12s" and re.search(r"large subunit ribosomal", definition, re.IGNORECASE):
        return {
            "resolved": {"12s"},  # the original name is kept, the conflict is recorded (as for 18s/28s)
            "reason": "12s matched with 'large subunit ribosomal' in definition, possible 16S mismatch, flagged",
            "conflict_reason": "12s conflicts with 'large subunit ribosomal' description",
        }
    if gene == "16s" and re.search(r"small subunit ribosomal", definition, re.IGNORECASE):
        return {
            "resolved": {"16s"},
            "reason": "16s matched with 'small subunit ribosomal' in definition, possible 12S mismatch, flagged",
            "conflict_reason": "16s conflicts with 'small subunit ribosomal' description",
        }
    if gene in ("18s", "28s") and "mitochondrial" in definition.lower():
        # Always preserve original gene name, record conflict in metadata
        if length_val > 800:
            return {
                "resolved": {gene},
                "reason": f"{gene} in definition with 'mitochondrial' but length {length_val}>800, kept {gene}",
                "conflict_reason": f"{gene} conflicts with mitochondrial annotation",
            }
        else:
            return {
                "resolved": {gene},
                "reason": f"{gene} in definition with 'mitochondrial' and length {length_val}<=800, "
                "flagged as uncertain",
                "conflict_reason": f"{gene} conflicts with mitochondrial annotation, short length",
            }
    return None


def conflict_check_and_reassign(
    gene_type_str: str,
    definition: str,
    topology: str,
    group: str,
    length_val: int,
    mt_min: int | None = 3000,
    cp_min: int | None = 20000,
) -> dict[str, Any]:
    genes = {g.strip().lower() for g in gene_type_str.split(",") if g.strip()}
    original = gene_type_str
    notes = []

    skip_types = {"mtgenome", "plastome", "18-28s", "its1-its2"}
    if genes & skip_types:
        notes.append("Special type (mtgenome/plastome/18-28s/its1-its2), skip conflict check")
        return _build_result(gene_type_str, "", notes, original=original)

    conflict = False
    conflict_reason = ""
    resolved = set(genes)

    if group.lower() == "mito":
        resolved, conflict, conflict_reason, notes = _conflict_mito(genes, definition, length_val, notes, mt_min)
    elif group.lower() == "plastid":
        resolved, conflict, conflict_reason, notes = _conflict_plastid(genes, length_val, notes, cp_min)
    else:
        resolved, conflict, conflict_reason, notes = _conflict_nuclear(genes, definition, length_val, notes)

    new_gene_type = ",".join(sorted(resolved))

    if len(resolved) > 1 and not conflict:
        conflict = True
        conflict_reason = f"Multiple candidate gene types: {', '.join(sorted(resolved))}"
        notes.append("Multiple candidates flagged as conflict")

    match_source = "reassign_mode" if new_gene_type != original else ""
    if new_gene_type != original:
        notes.append(f"Reassigned: [{original}] -> [{new_gene_type}]")

    return _build_result(
        new_gene_type, match_source, notes, original=original, conflict=conflict, conflict_reason=conflict_reason
    )


def _conflict_mito(
    genes: set[str], definition: str, length_val: int, notes: list[str], mt_min: int | None = 3000
) -> tuple[set[str], bool, str, list[str]]:
    conflict = False
    conflict_reason = ""
    resolved = set(genes)

    if len(genes) > 2:
        # a long multi-gene record is a mitogenome only if its title says so (plant mito fragments of 5-30 kb
        # list several genes: Ginkgo MW354233); otherwise the genes are kept and flagged
        if (mt_min is None or length_val >= mt_min) and "genome" in definition.lower():
            resolved = {"mtgenome"}
            notes.append(f"Mitochondrial multi-match ({len(genes)}) + Length>={mt_min} -> reassigned to mtgenome")
        else:
            conflict = True
            conflict_reason = f"Mitochondrial multi-match ({len(genes)}) but length {length_val}<{mt_min}"
            notes.append("Mitochondrial multi-match but insufficient length, kept original + flagged")
        return resolved, conflict, conflict_reason, notes

    if len(genes) == 2:
        has_trna = re.search(r"trna", definition, re.IGNORECASE)
        has_complete = re.search(r"complete sequence", definition, re.IGNORECASE)
        if has_trna or has_complete:
            conflict = True
            trigger = "tRNA" if has_trna else "complete sequence"
            conflict_reason = f"Two mitochondrial genes + definition contains '{trigger}'"
            notes.append(f"Dual match + '{trigger}' in definition, kept original + flagged")
            return resolved, conflict, conflict_reason, notes

    for gene in sorted(genes):
        result = _check_ribosomal_conflict(gene, definition, length_val)
        if result:
            conflict = True
            conflict_reason = result["conflict_reason"]
            notes.append(result["reason"])
            return resolved, conflict, conflict_reason, notes

    return resolved, conflict, conflict_reason, notes


def _conflict_plastid(
    genes: set[str], length_val: int, notes: list[str], cp_min: int | None = 20000
) -> tuple[set[str], bool, str, list[str]]:
    """Like _conflict_mito: more than two plastid genes in a record of plastome size -> plastome."""
    if len(genes) > 2 and (cp_min is None or length_val >= cp_min):
        notes.append(f"Plastid multi-match ({len(genes)}) + Length>={cp_min} -> reassigned to plastome")
        return {"plastome"}, False, "", notes
    return set(genes), False, "", notes


def _conflict_nuclear(
    genes: set[str], definition: str, length_val: int, notes: list[str]
) -> tuple[set[str], bool, str, list[str]]:
    conflict = False
    conflict_reason = ""
    resolved = set(genes)

    for gene in sorted(genes):
        result = _check_ribosomal_conflict(gene, definition, length_val)
        if result:
            conflict = True
            conflict_reason = result["conflict_reason"]
            notes.append(result["reason"])
            return resolved, conflict, conflict_reason, notes

    return resolved, conflict, conflict_reason, notes


# =========================
# Round 2: Recheck unmatched
# =========================
def recheck_match_row(row: dict[str, Any], config: MatchConfig) -> dict[str, Any]:
    definition = clean_definition(row.get("Definition", ""))
    organelle = str(row.get("organelle", "")).strip()
    molecule_type = str(row.get("MoleculeType", "")).lower().strip()
    topology = str(row.get("Topology", "")).lower().strip()
    length_val = clean_length(row.get("Length", "0"))
    pcr_primers = str(row.get("PCR_primers", "")).strip() if "PCR_primers" in row else ""

    candidates = set()
    notes = []
    conflict_notes = []

    mito_terms = ["mitochondrion", "mitochondrial", "mitochondria"]
    has_mito_def = any(t in definition for t in mito_terms)
    has_mito_org = any(t in organelle.lower() for t in mito_terms)

    group = _record_group(organelle.lower(), definition)
    if topology == "circular":
        genome = "plastome" if group == "plastid" else "mtgenome"
        candidates.add(genome)
        notes.append(f"Topology=circular ({group}) -> direct {genome}")
    else:
        for gene, keywords in config.gene_dict.items():
            if not keywords:
                continue
            if not any_match(keywords, definition):
                continue
            if _negative_hit(gene, definition):
                conflict_notes.append(f"{gene}: keyword matched but excluded by context")
                continue

            if gene == "mtgenome":
                mt_min = _mtgenome_min(config)
                if (
                    (mt_min is None or length_val >= mt_min)
                    and molecule_type == "dna"
                    and (has_mito_def or has_mito_org)
                ):
                    candidates.add("mtgenome")
                    notes.append(f"mtgenome: keyword + Length>={mt_min} + DNA + mitochondrial features")
                else:
                    conflict_notes.append(
                        f"mtgenome: conditions not met (Length={length_val}, "
                        f"MoleculeType={molecule_type}, "
                        f"mito_in_def={has_mito_def}, mito_in_org={has_mito_org})"
                    )
            else:
                candidates.add(gene)
                notes.append(f"{gene}: keyword matched")

        if not candidates and pcr_primers:
            coi_primers = ["HCO2198", "lc01490", "LCO1490"]
            if any(re.search(rf"(?:fwd|rev)_name:\s*{re.escape(p)}", pcr_primers, re.IGNORECASE) for p in coi_primers):
                candidates.add("coi")
                notes.append("PCR primer match (HCO2198/lc01490) -> coi")

            if re.search(r"(?:fwd|rev)_name:\s*28SC1", pcr_primers, re.IGNORECASE):
                candidates.add("28s")
                notes.append("PCR primer match (28SC1) -> 28s")

    # plants have rRNA genes in all three compartments (mitochondrial rrn26 is not nuclear 28S, plastid 16S is not
    # mitochondrial 16S): there a gene of another compartment is dropped. In animals the record keeps the gene and
    # is flagged, because "18S ...; mitochondrial" is a submitter's slip, not a mitochondrial 18S.
    plant = group == "plastid" or "viridiplantae" in str(row.get("Taxonomy", "")).lower()
    other = (
        (config.nuclear_genes | config.plastid_genes) - {"mtgenome"}
        if plant and group == "mito"
        else (config.mito_genes | config.nuclear_genes) - {"plastome"}
        if group == "plastid"
        else set()
    )
    dropped = sorted(candidates & other)
    if dropped:
        candidates -= set(dropped)
        notes.append(f"dropped {', '.join(dropped)}: gene of another compartment than this {group} record")

    if candidates:
        if len(candidates) > 1:
            conflict_notes.append(f"Multiple candidates: {', '.join(sorted(candidates))}")
        elif len(candidates) == 1:
            gene_assigned = next(iter(candidates))
            if gene_assigned in config.mito_genes and not (has_mito_def or has_mito_org):
                conflict_notes.append(
                    f"Mitochondrial candidate ({gene_assigned}) but no mitochondrial marker in organelle/definition"
                )
            if gene_assigned in config.nuclear_genes and (has_mito_def or has_mito_org):
                conflict_notes.append(
                    f"Nuclear candidate ({gene_assigned}) but mitochondrial marker in organelle/definition"
                )

    gene_type = ",".join(sorted(candidates))
    match_source = "recheck_mode" if gene_type else ""
    conflict = bool(conflict_notes) and bool(gene_type)

    return _build_result(
        gene_type,
        match_source,
        notes,
        original=gene_type,
        conflict=conflict,
        conflict_reason="; ".join(conflict_notes) if conflict else "",
    )


# =========================
# Result builder
# =========================
def _build_result(
    gene_type: str,
    match_source: str,
    notes: list[str],
    original: str = "",
    conflict: bool = False,
    conflict_reason: str = "",
) -> dict[str, Any]:
    return {
        "gene_type": gene_type,
        "match_source": match_source,
        "Original_match": original or gene_type,
        "conflict": conflict,
        "Conflict_reason": conflict_reason,
        "Assignment_reason": "; ".join(notes),
    }


# =========================
# Row processing
# =========================
def process_row(row: dict[str, Any], config: MatchConfig) -> dict[str, Any]:
    prelim = first_round_match(row, config)
    gene_type_str = prelim["gene_type"]

    if not gene_type_str:
        return prelim

    definition = clean_definition(row.get("Definition", ""))
    group = _record_group(str(row.get("organelle", "")).lower(), definition)
    length_val = clean_length(row.get("Length", "0"))

    final = conflict_check_and_reassign(
        gene_type_str,
        definition,
        "",
        group,
        length_val,
        mt_min=_mtgenome_min(config),
        cp_min=parse_interval(config.get_gene_length_range("plastome"))[0],
    )

    if final["match_source"] == "reassign_mode":
        match_source = "first_match -> conflict_check"
    else:
        match_source = prelim["match_source"]

    reasons = []
    if prelim["Assignment_reason"]:
        reasons.append(prelim["Assignment_reason"])
    if final["Assignment_reason"]:
        reasons.append(f"[conflict check] {final['Assignment_reason']}")

    return {
        "gene_type": final["gene_type"],
        "match_source": match_source,
        "Original_match": final["Original_match"],
        "conflict": final["conflict"],
        "Conflict_reason": final["Conflict_reason"],
        "Assignment_reason": "; ".join(reasons),
    }


# =========================
# Main flows
# =========================
def process_prematch(
    input_file: str,
    output_dir: str,
    config: MatchConfig,
    prematch_output: str,
    multip_output: str,
    unmatched_output: str,
    conflicted_output: str,
    moleculetype: str,
    mol_type: str,
    length_filter: str,
    organelle_filter: str,
) -> pd.DataFrame:
    os.makedirs(output_dir, exist_ok=True)

    ext = os.path.splitext(input_file)[1].lower()
    if ext not in (".csv", ".tsv", ".xls", ".xlsx", ".txt"):
        logger.error(f"Unsupported format: {ext}")
        raise UnsupportedFormatError(f"Unsupported format: {ext}")
    df = utils_read_table(input_file)
    df.columns = df.columns.str.strip()

    df["g2t_length_raw"] = df["Length"].astype(str)
    df["Length"] = df["Length"].astype(str).str.replace(r"\s*bp\s*", "", regex=True).str.replace(",", "")
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
    df["g2t_input_row"] = range(len(df))
    removed: list[dict] = []

    missing = df["Length"].isna()
    removed.extend(_filtered_rows(df[missing], "classify round 1", "length missing or not numeric"))
    df = df[~missing]
    lo, hi = parse_interval(config.global_length_range)
    keep = pd.Series(True, index=df.index)
    if lo is not None:
        keep &= df["Length"] >= lo
    if hi is not None:
        keep &= df["Length"] <= hi
    # plastid records (plastomes ~120-220 kb) are checked against --length2_cpgenes instead
    plastid = pd.Series(
        [
            _record_group(str(o).lower() if pd.notna(o) else "", str(d) if pd.notna(d) else "") == "plastid"
            for o, d in zip(
                df.get("organelle", pd.Series("", index=df.index)), df.get("Definition", pd.Series("", index=df.index))
            )
        ],
        index=df.index,
        dtype=bool,
    )
    cp_lo, cp_hi = parse_interval(config.length2_cpgenes)
    cp_ok = plastid & (df["Length"] >= (cp_lo or 0)) & (df["Length"] <= (cp_hi if cp_hi is not None else float("inf")))
    keep |= cp_ok
    removed.extend(
        _filtered_rows(
            df[~keep],
            "classify round 1",
            lambda r: (
                f"length {int(r['Length'])} bp outside {config.global_length_range} "
                "(global length filter; change with --length_range2_all)"
            ),
        )
    )
    df = df[keep]
    logger.info(f"After global length filter: {len(df)} records")

    df = _filter_dataframe_logged(df, moleculetype, mol_type, length_filter, organelle_filter, removed)
    logger.info(f"After additional filters: {len(df)} records")

    results = []
    conflict_records = []
    for row_dict in df.to_dict("records"):  # itertuples renamed columns with spaces
        res = process_row(row_dict, config)
        row_dict["gene_type"] = res["gene_type"]
        row_dict["match_source"] = res["match_source"]
        row_dict["Original_match"] = res["Original_match"]
        row_dict["Conflict"] = res["conflict"]
        row_dict["Conflict_reason"] = res["Conflict_reason"]
        row_dict["Assignment_reason"] = res["Assignment_reason"]
        if res["conflict"]:
            conflict_records.append(row_dict)
        results.append(row_dict)

    res_df = pd.DataFrame(results) if results else pd.DataFrame(columns=list(df.columns) + DIAG_COLUMNS)
    if "LocusID" in res_df.columns and len(res_df):
        dup = res_df.duplicated(subset=["LocusID"], keep="first")
        removed.extend(_filtered_rows(res_df[dup], "classify round 1", "duplicate LocusID (first occurrence kept)"))
        res_df = res_df[~dup]
    _write_filtered(output_dir, removed, append=False)
    if removed:
        logger.info(f"Filtered out {len(removed)} records -> {FILTERED_FILE} (with reasons)")
    res_df = res_df.drop(columns=["g2t_input_row", "g2t_length_raw"], errors="ignore")

    if not config.add_assignment_reasons and "Assignment_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Assignment_reason"])
    if not config.record_conflict_reasons and "Conflict_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Conflict_reason"])

    res_df = reorder_columns(res_df)

    assigned_df = res_df[res_df["gene_type"] != ""]
    unmatched_df = res_df[res_df["gene_type"] == ""]

    _save_csv(assigned_df, output_dir, prematch_output, "Matched")
    _save_csv(unmatched_df, output_dir, unmatched_output, "Unmatched")

    if conflict_records:
        conflict_df = reorder_columns(
            pd.DataFrame(conflict_records).drop(columns=["g2t_input_row", "g2t_length_raw"], errors="ignore")
        )
        _save_csv(conflict_df, output_dir, conflicted_output, "Conflict")

    multip_df = assigned_df[assigned_df["gene_type"].apply(lambda x: len(x.split(",")) > 1 if x else False)]
    if not multip_df.empty:
        _save_csv(multip_df, output_dir, multip_output, "Multi-match")

    total = len(res_df)
    logger.info(
        f"Round 1 summary: total={total}, matched={len(assigned_df)}, "
        f"conflict={len(conflict_records)}, unmatched={len(unmatched_df)}"
    )
    return res_df


def process_recheck(
    input_file: str,
    output_dir: str,
    config: MatchConfig,
    apply_length_filter: bool = True,
    assigned_out: str = "assigned_genes_types2.csv",
    multip_out: str = "multiple_matched2.csv",
    unmatched_out: str = "unmatched_sequences2.csv",
) -> pd.DataFrame:
    os.makedirs(output_dir, exist_ok=True)

    ext = os.path.splitext(input_file)[1].lower()
    if ext not in (".csv", ".tsv", ".xls", ".xlsx", ".txt"):
        logger.error(f"Unsupported format: {ext}")
        raise UnsupportedFormatError(f"Unsupported format: {ext}")
    df = utils_read_table(input_file)
    df.columns = df.columns.str.strip()

    if df.empty:
        logger.info("Round 2: no unmatched records to recheck")
        _save_csv(df, output_dir, assigned_out, "Round 2 matched")
        _save_csv(df, output_dir, unmatched_out, "Round 2 unmatched")
        return df

    # Use union of mt+nt ranges instead of nuclear-only
    if apply_length_filter:
        df["Length"] = df["Length"].astype(str).str.replace(r"\s*bp\s*", "", regex=True).str.replace(",", "")
        df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
        ranges = [parse_interval(r) for r in (config.length2_mtgenes, config.length2_ntgenes, config.length2_cpgenes)]
        lo = min(r[0] or 0 for r in ranges)
        hi = max(r[1] or float("inf") for r in ranges)
        if lo is not None:
            df = df[df["Length"] >= lo]
        if hi is not None and hi != float("inf"):
            df = df[df["Length"] <= hi]

    results = []
    conflict_records = []
    for row_dict in df.to_dict("records"):  # itertuples renamed columns with spaces
        res = recheck_match_row(row_dict, config)
        row_dict["gene_type"] = res["gene_type"]
        row_dict["match_source"] = res["match_source"]
        row_dict["Original_match"] = res["Original_match"]
        row_dict["Conflict"] = res["conflict"]
        row_dict["Conflict_reason"] = res["Conflict_reason"]
        row_dict["Assignment_reason"] = res["Assignment_reason"]
        if res["conflict"]:
            conflict_records.append(row_dict)
        results.append(row_dict)

    res_df = pd.DataFrame(results) if results else pd.DataFrame(columns=list(df.columns) + DIAG_COLUMNS)
    if "LocusID" in res_df.columns:
        res_df = res_df.drop_duplicates(subset=["LocusID"])

    if not config.add_assignment_reasons and "Assignment_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Assignment_reason"])
    if not config.record_conflict_reasons and "Conflict_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Conflict_reason"])

    res_df = reorder_columns(res_df)

    assigned_df = res_df[res_df["gene_type"] != ""]
    unmatched_df = res_df[res_df["gene_type"] == ""]
    multi_df = assigned_df[assigned_df["gene_type"].apply(lambda x: len(x.split(",")) > 1 if x else False)]

    _save_csv(assigned_df, output_dir, assigned_out, "Round 2 matched")
    _save_csv(unmatched_df, output_dir, unmatched_out, "Round 2 unmatched")
    if not multi_df.empty:
        _save_csv(multi_df, output_dir, multip_out, "Round 2 multi-match")

    if conflict_records:
        conflict_file = os.path.join(output_dir, "conflicted_genes.csv")
        conflict_df = reorder_columns(pd.DataFrame(conflict_records))
        if os.path.exists(conflict_file):
            try:
                existing = pd.read_csv(conflict_file, dtype=str)
                combined = pd.concat([existing, conflict_df]).drop_duplicates()
                combined.to_csv(conflict_file, index=False)
                logger.info(f"Appended {len(conflict_df)} conflict records to {conflict_file}")
            except Exception as e:
                logger.error(f"Failed to merge conflict records: {e}")
        else:
            _save_csv(conflict_df, output_dir, "conflicted_genes.csv", "Round 2 conflict")

    logger.info(
        f"Round 2 summary: total={len(res_df)}, matched={len(assigned_df)}, "
        f"conflict={len(conflict_records)}, unmatched={len(unmatched_df)}"
    )
    return res_df


def extract_genes(df: pd.DataFrame, gene_types: list[str] | None, output_dir: str) -> None:
    if not gene_types:
        all_genes: set[str] = set()
        for gt in df["gene_type"].dropna():
            all_genes.update(s.strip() for s in str(gt).split(",") if s.strip())
        gene_types = sorted(all_genes)
        logger.info(f"Auto-detected gene types: {gene_types}")

    for gene in gene_types:
        subset = df[
            df["gene_type"].apply(lambda x, g=gene: g in [s.strip() for s in str(x).split(",")] if x else False)
        ]
        if subset.empty:
            logger.warning(f"No records for gene type '{gene}'")
            continue
        safe_name = re.sub(r"[^\w\-]", "_", gene)
        out_file = os.path.join(output_dir, f"{safe_name}_extracted.csv")
        try:
            reorder_columns(subset).to_csv(out_file, index=False)
            logger.info(f"  {gene} -> {out_file} ({len(subset)} records)")
        except Exception as e:
            logger.error(f"Failed to save {gene} extraction: {e}")


def _save_csv(df: pd.DataFrame, output_dir: str, filename: str, label: str) -> None:
    path = os.path.join(output_dir, filename)
    try:
        df.to_csv(path, index=False)
        logger.info(f"{label}: {path} ({len(df)} records)")
    except Exception as e:
        logger.error(f"Failed to save {label}: {e}")


def categorize_unmatched_definition(def_text: str) -> str:
    """Categorize unmatched sequence by its Definition field."""
    def_lower = str(def_text).lower()

    # Genome/assembly data
    if "whole genome shotgun" in def_lower:
        return "whole_genome_shotgun"
    elif "genome assembly" in def_lower:
        return "genome_assembly"
    elif "tpa_asm" in def_lower or "tpa:" in def_lower:
        return "TPA_assembly"
    elif "chromosome:" in def_lower:
        return "chromosome_data"
    elif "scaffold" in def_lower:
        return "scaffold"
    elif "contig" in def_lower:
        return "contig"

    # RNA types
    elif "mrna sequence" in def_lower or "cdna clone" in def_lower:
        return "mRNA/cDNA"
    elif "trna" in def_lower:
        return "tRNA"
    elif "rrna" in def_lower and "mitochondrial" not in def_lower:
        return "nuclear_rRNA"

    # Other genes not in standard list
    elif "nad" in def_lower or "nadh" in def_lower:
        return "NADH_dehydrogenase"
    elif "atp synthase" in def_lower or "atpase" in def_lower:
        return "ATP_synthase"
    elif "creatine kinase" in def_lower:
        return "creatine_kinase"
    elif "cytochrome b" in def_lower or "cytb" in def_lower:
        return "cytochrome_b_variant"
    else:
        return "other"


def generate_unmatched_report(unmatched_df: pd.DataFrame, output_dir: str) -> str:
    """Generate a summary report of unmatched sequences by category."""
    if unmatched_df.empty:
        return ""

    # Categorize each unmatched sequence
    unmatched_df = unmatched_df.copy()
    unmatched_df["category"] = unmatched_df["Definition"].apply(categorize_unmatched_definition)

    # Generate summary
    summary_lines = []
    summary_lines.append("=" * 60)
    summary_lines.append("UNMATCHED SEQUENCES REPORT")
    summary_lines.append("=" * 60)
    summary_lines.append(f"Total unmatched: {len(unmatched_df)}")
    summary_lines.append("")

    category_counts = unmatched_df["category"].value_counts()
    summary_lines.append("Category breakdown:")
    summary_lines.append("-" * 40)

    for cat, count in category_counts.items():
        pct = count / len(unmatched_df) * 100
        summary_lines.append(f"  {cat}: {count} ({pct:.1f}%)")

    summary_lines.append("")
    summary_lines.append("Sample definitions per category:")
    summary_lines.append("-" * 40)

    for cat in category_counts.index[:5]:  # Top 5 categories
        samples = unmatched_df[unmatched_df["category"] == cat]["Definition"].head(3).tolist()
        summary_lines.append(f"\n[{cat}]")
        for s in samples:
            summary_lines.append(f"  - {str(s)[:100]}...")

    summary_lines.append("")
    summary_lines.append("=" * 60)

    # Save report
    report_path = os.path.join(output_dir, "unmatched_report.txt")
    report_content = "\n".join(summary_lines)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_content)

    # Also save categorized CSV
    categorized_path = os.path.join(output_dir, "unmatched_categorized.csv")
    unmatched_df.to_csv(categorized_path, index=False)

    logger.info(f"Unmatched report saved: {report_path}")
    logger.info(f"Categorized unmatched: {categorized_path}")

    # Print summary to console
    print("\n" + report_content)

    return report_path


# =========================
# Main
# =========================
def _run_cli(args: argparse.Namespace) -> None:
    os.makedirs(args.output, exist_ok=True)
    pkg_logger = logging.getLogger("g2t")
    old_level = pkg_logger.level
    handler = setup_logging(os.path.join(args.output, "classify.log"))
    try:
        _run_cli_steps(args)
    finally:
        pkg_logger.removeHandler(handler)
        pkg_logger.setLevel(old_level)
        handler.close()


def _run_cli_steps(args: argparse.Namespace) -> None:
    logger.info("=" * 60)
    logger.info("Gene Type Classifier (refactored)")
    logger.info("=" * 60)

    output_dir = args.output

    mtgenes = [s.strip() for s in args.mtgenes_list.split(",") if s.strip()]
    ntgenes = [s.strip() for s in args.ntgenes_list.split(",") if s.strip()]
    cpgenes = [s.strip() for s in args.cpgenes_list.split(",") if s.strip()]
    gene_dict = load_path_dict(args.path_dict) if args.path_dict else default_gene_dict()
    mito_genes, nuclear_genes, plastid_genes = init_gene_groups(mtgenes, ntgenes, cpgenes, gene_dict)

    gene_ranges = {
        "mtgenome": args.length_range2_mtgenome,
        "plastome": args.length_range2_plastome,
        "18-28s": args.length_range2_18_28s,
        "coi": args.length_range2_coi,
        "its1-its2": args.length_range2_its,
        "cob": args.length_range2_cob,
        "12s": args.length_range2_12s,
        "16s": args.length_range2_16s,
        "18s": args.length_range2_18s,
        "28s": args.length_range2_28s,
        "h3": args.length_range2_h3,
    }
    if args.gene_length_range:
        for item in args.gene_length_range:
            if "=" in item:
                gene, rng = item.split("=", 1)
                gene_ranges[gene.strip().lower()] = rng.strip()

    config = MatchConfig(
        if_remove_matched=args.if_remove_matched_sequences,
        length2_mtgenes=args.length2_mtgenes,
        length2_ntgenes=args.length2_ntgenes,
        length2_cpgenes=args.length2_cpgenes,
        global_length_range=args.length_range2_all,
        gene_length_ranges=gene_ranges,
        add_assignment_reasons=args.if_add_assigned_gene_types_reasons.lower() == "yes",
        record_conflict_reasons=args.if_recorded_conflict_reasons.lower() == "yes",
        mito_genes=mito_genes,
        nuclear_genes=nuclear_genes,
        plastid_genes=plastid_genes,
        gene_dict=gene_dict,
    )

    gene_types_arg = args.which_gene_types_extract
    result = classify(
        args.input,
        output_dir,
        config,
        prematch_output=args.prematch_output_file_name,
        unmatched_output=args.unmatched_sequences,
        if_recheck=args.if_recheck_unmatched.lower() == "yes",
        apply_length_filter_recheck=args.if_length_range_work2_recheck.lower() == "yes",
        extract_gene_types=[s.strip().lower() for s in gene_types_arg.split(",")] if gene_types_arg else None,
        multiple_output=args.multip_matched_genestype,
        conflicted_output=args.conflicted_sequences,
        moleculetype=args.moleculetype,
        mol_type=args.mol_type,
        length=args.length,
        organelle=args.organelle,
    )
    if not result.success:
        sys.exit(1)
    logger.info("=" * 60)
    logger.info("Processing complete")
    logger.info("=" * 60)


def build_record_status(
    input_file: str,
    output_dir: str,
    if_recheck: bool,
    config: MatchConfig,
    unmatched_output: str = "unmatched_sequences.csv",
    prematch_output: str = "assigned_genes_type.csv",
) -> pd.DataFrame:
    """One row per input row: assigned (gene_type) / unmatched / filtered, each with a reason."""
    src = utils_read_table(input_file)
    src.columns = src.columns.str.strip()

    def read(name):
        path = os.path.join(output_dir, name)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            try:
                return pd.read_csv(path, dtype=str).fillna("")
            except pd.errors.EmptyDataError:
                pass
        return pd.DataFrame()

    filt = read(FILTERED_FILE)
    filt_reason = {str(r["input_row"]): r["reason"] for _, r in filt.iterrows()} if len(filt) else {}
    assigned = read("assigned_genes_types_all.csv")
    if assigned.empty:
        assigned = read(prematch_output)
    gene_of = dict(zip(assigned.get("LocusID", []), assigned.get("gene_type", [])))
    reason_of = dict(zip(assigned.get("LocusID", []), assigned.get("Assignment_reason", [""] * len(assigned))))
    recheck_ran = if_recheck and os.path.exists(os.path.join(output_dir, "unmatched_sequences2.csv"))
    un2 = read("unmatched_sequences2.csv") if recheck_ran else pd.DataFrame()
    un2_ids = set(un2.get("LocusID", []))
    rows = []
    for i, r in src.iterrows():
        lid = str(r.get("LocusID", ""))
        base = {
            "input_row": i,
            "ACCESSION": r.get("ACCESSION", ""),
            "LocusID": lid,
            "Organism": r.get("Organism", r.get("organism", "")),
            "Length": r.get("Length", ""),
            "Definition": r.get("Definition", ""),
        }
        if str(i) in filt_reason:
            rows.append({**base, "status": "filtered", "gene_type": "", "reason": filt_reason[str(i)]})
        elif lid in gene_of:
            rows.append({**base, "status": "assigned", "gene_type": gene_of[lid], "reason": reason_of.get(lid, "")})
        else:
            reason = "no gene type matched (DEFINITION keywords, round 1"
            if recheck_ran and lid not in un2_ids:
                reason += (
                    f"); skipped in round-2 recheck: length outside {config.length2_mtgenes} / {config.length2_ntgenes}"
                    f" / {config.length2_cpgenes}"
                )
            else:
                reason += " and 2)" if recheck_ran else ")"
            rows.append({**base, "status": "unmatched", "gene_type": "", "reason": reason})
    out = pd.DataFrame(
        rows,
        columns=[
            "input_row",
            "ACCESSION",
            "LocusID",
            "Organism",
            "Length",
            "Definition",
            "status",
            "gene_type",
            "reason",
        ],
    )
    out.to_csv(os.path.join(output_dir, STATUS_FILE), index=False)
    counts = out["status"].value_counts().to_dict()
    logger.info(
        f"Record status: input={len(out)}, assigned={counts.get('assigned', 0)}, "
        f"unmatched={counts.get('unmatched', 0)}, filtered={counts.get('filtered', 0)} -> {STATUS_FILE}"
    )
    return out


def classify(
    input_file: str,
    output_dir: str,
    config: MatchConfig | None = None,
    prematch_output: str = "assigned_genes_type.csv",
    unmatched_output: str = "unmatched_sequences.csv",
    if_recheck: bool = True,
    apply_length_filter_recheck: bool = True,
    extract_gene_types: list | None = None,
    *,
    multiple_output: str = "multiple_matched.csv",
    conflicted_output: str = "conflicted_genes.csv",
    moleculetype: str = "all",
    mol_type: str = "all",
    length: str = "all",
    organelle: str = "all",
) -> StepResult:
    """Classify gene types from metadata CSV. Returns StepResult. ``g2t-classify`` calls this too.

    The keyword-only arguments are the file names and record filters of ``g2t-classify``.
    """
    import time

    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    if config is None:
        config = MatchConfig()
    _clean_outputs(output_dir, [prematch_output, unmatched_output, multiple_output, conflicted_output])

    try:
        logger.info("-" * 40)
        logger.info("Round 1: Gene type matching")
        logger.info("-" * 40)
        prematch_df = process_prematch(
            input_file,
            output_dir,
            config,
            prematch_output,
            multiple_output,
            unmatched_output,
            conflicted_output,
            moleculetype,
            mol_type,
            length,
            organelle,
        )
        if extract_gene_types:
            extract_genes(prematch_df, extract_gene_types, output_dir)

        unmatched_file = os.path.join(output_dir, unmatched_output)
        if if_recheck:
            if os.path.exists(unmatched_file) and os.path.getsize(unmatched_file) > 50:
                logger.info("-" * 40)
                logger.info("Round 2: Unmatched record recheck")
                logger.info("-" * 40)
                process_recheck(
                    unmatched_file,
                    output_dir,
                    config,
                    apply_length_filter=apply_length_filter_recheck,
                )
            else:
                logger.warning(f"Unmatched file not found or empty: {unmatched_file}, skipping round 2")
        all_file = _write_all(output_dir, prematch_output)
        logger.info(f"Merged results: {all_file}")

        try:
            build_record_status(input_file, output_dir, if_recheck, config, unmatched_output, prematch_output)
        except Exception as e:  # noqa: BLE001 - accounting must not break classification
            logger.warning(f"Failed to write {STATUS_FILE}: {e}")

        if os.path.exists(unmatched_file):
            try:
                unmatched_df = pd.read_csv(unmatched_file, dtype=str)
                if not unmatched_df.empty:
                    generate_unmatched_report(unmatched_df, output_dir)
            except (OSError, pd.errors.EmptyDataError) as e:
                logger.warning(f"Failed to generate unmatched report: {e}")

        rows = len(pd.read_csv(all_file, dtype=str))
        return StepResult(
            success=True,
            output_file=all_file,
            rows=rows,
            elapsed=time.time() - start_time,
        )
    except Exception as e:
        logger.exception(f"Classify failed: {e}")
        return StepResult(success=False, output_file="", elapsed=time.time() - start_time)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Gene Type Classifier: GenBank metadata -> gene_type labels",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-i", "--input", required=True, help="Input file (CSV/TSV/Excel)")
    parser.add_argument("-o", "--output", default="assign_genes_type", help="Output directory")

    parser.add_argument(
        "--prematch_output_file_name", default="assigned_genes_type.csv", help="Round-1 assignments (file name in -o)"
    )
    parser.add_argument(
        "--multip_matched_genestype", default="multiple_matched.csv", help="Records matching several gene types"
    )
    parser.add_argument(
        "--unmatched_sequences", default="unmatched_sequences.csv", help="Records matching no gene type"
    )
    parser.add_argument(
        "--conflicted_sequences", default="conflicted_genes.csv", help="Records whose title contradicts the gene type"
    )

    parser.add_argument(
        "--if_remove_matched_sequences",
        "-ifrms",
        default="no",
        help="yes: stop at the first matching gene type (one gene per record)",
    )

    parser.add_argument("--length_range2_mtgenome", default="3000:500000", help="Length window (bp) of mitogenomes")
    parser.add_argument("--length_range2_plastome", default="20000:300000", help="Length window (bp) of plastomes")
    parser.add_argument(
        "--length_range2_18_28s", default="5000:15000", help="Length window (bp) of whole 18S-ITS-28S records"
    )
    parser.add_argument("--length_range2_coi", default="none", help="Length window (bp) of coi; none = no limit")
    parser.add_argument("--length_range2_its", default="none", help="Length window (bp) of its; none = no limit")
    parser.add_argument("--length_range2_cob", default="none", help="Length window (bp) of cob; none = no limit")
    parser.add_argument("--length_range2_12s", default="none", help="Length window (bp) of 12s; none = no limit")
    parser.add_argument("--length_range2_16s", default="none", help="Length window (bp) of 16s; none = no limit")
    parser.add_argument("--length_range2_18s", default="none", help="Length window (bp) of 18s; none = no limit")
    parser.add_argument("--length_range2_28s", default="none", help="Length window (bp) of 28s; none = no limit")
    parser.add_argument("--length_range2_h3", default="none", help="Length window (bp) of h3; none = no limit")

    parser.add_argument("--gene_length_range", nargs="*", default=None, help="Per-gene length ranges: gene=lower:upper")

    parser.add_argument(
        "--length2_mtgenes", default="150:50000", help="Length (bp) a record needs to be typed as a mitochondrial gene"
    )
    parser.add_argument(
        "--length2_ntgenes", default="150:50000", help="Length (bp) a record needs to be typed as a nuclear gene"
    )
    parser.add_argument(
        "--length2_cpgenes", default="150:300000", help="Length (bp) a record needs to be typed as a plastid gene"
    )
    parser.add_argument(
        "--length_range2_all",
        default="150:50000",
        help="Global length filter of round 1 (plastid records use --length2_cpgenes); others are listed as filtered",
    )
    parser.add_argument("--length", default="all", help="Extra length filter, e.g. 300:2000 (all = none)")
    parser.add_argument("--organelle", default="all", help="Keep only this /organelle value, e.g. mitochondrion")
    parser.add_argument("--moleculetype", default="all", help="Keep only this LOCUS molecule type, e.g. DNA")
    parser.add_argument("--mol_type", default="all", help="Keep only these /mol_type values (comma-separated)")

    parser.add_argument(
        "--mtgenes_list",
        default=",".join(MITO_TYPES),
        help="Mitochondrial gene types searched",
    )
    parser.add_argument("--ntgenes_list", default=",".join(NUCLEAR_TYPES), help="Nuclear gene types searched")
    parser.add_argument(
        "--cpgenes_list",
        default=",".join(PLASTID_TYPES),
        help="Plastid gene types searched",
    )
    parser.add_argument("--path_dict", default=None, help="External gene dictionary file")

    # "--which_gene_types_extact" (misspelt) is still accepted
    parser.add_argument(
        "--which_gene_types_extract",
        "--which_gene_types_extact",
        "-wgte",
        default=None,
        help="Also write one file per gene type for these types (comma-separated)",
    )

    parser.add_argument(
        "--if_recheck_unmatched",
        default="yes",
        help="yes: round 2 re-checks unmatched records against all compartments (circular = genome)",
    )
    parser.add_argument(
        "--if_length_range_work2_recheck", default="yes", help="yes: round 2 applies the length windows too"
    )
    parser.add_argument(
        "--if_add_assigned_gene_types_reasons", default="yes", help="yes: write the Assignment_reason column"
    )
    parser.add_argument("--if_recorded_conflict_reasons", default="yes", help="yes: write the Conflict_reason column")

    args = parser.parse_args(argv)

    check_argv = argv if argv else sys.argv
    if "--which_gene_types_extact" in check_argv:
        warnings.warn(
            "--which_gene_types_extact is deprecated (typo), use --which_gene_types_extract",
            DeprecationWarning,
            stacklevel=2,
        )

    _run_cli(args)


if __name__ == "__main__":
    main()
