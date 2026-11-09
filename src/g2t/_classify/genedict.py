"""The gene dictionary: packaged synonyms, user dictionaries, keyword patterns."""

from __future__ import annotations

import importlib.util
import logging
import os
import re
from functools import cache, lru_cache

import pandas as pd

from g2t.genes import MITO_TYPES, NUCLEAR_TYPES, PLASTID_TYPES

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


# =========================
# Gene dictionary: the packaged data/gene_dict.yaml is the single source of synonyms
# =========================
def _yaml_path() -> str:
    try:
        from importlib.resources import files

        return str(files("g2t") / "data" / "gene_dict.yaml")
    except Exception:  # noqa: BLE001 - zipped or unusual installs
        return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "gene_dict.yaml")


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
