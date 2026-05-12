import re
import pandas as pd
import argparse
import logging
import os
import sys
import warnings
import importlib.util
from dataclasses import dataclass, field
from typing import Tuple, Optional, Dict, Any, List, Set
from g2t.utils import StepResult, parse_interval, length_in_range

# =========================
# Logging configuration
# =========================
class CustomFormatter(logging.Formatter):
    def format(self, record):
        if record.levelno == logging.DEBUG:
            record.levelname = "CHECKING"
        return super().format(record)


def setup_logging(log_file: str = 'extract_methods.log') -> None:
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = CustomFormatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


setup_logging()
logger = logging.getLogger(__name__)


# =========================
# GENE_DICT: load from YAML or use built-in fallback
# =========================
def _load_gene_dict_yaml(yaml_path=None):
    """Load gene dictionary from YAML file. Returns None if unavailable."""
    if yaml_path is None:
        try:
            from importlib.resources import files
            yaml_path = str(files("g2t") / "data" / "gene_dict.yaml")
        except Exception:
            yaml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gene_dict.yaml")
    if not os.path.exists(yaml_path):
        return None
    try:
        import yaml
        with open(yaml_path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
        result = {}
        for k, v in data.get('gene_types', {}).items():
            synonyms = v.get('synonyms', []) if isinstance(v, dict) else (v if isinstance(v, list) else [])
            result[k.lower()] = list({str(s).lower() for s in synonyms if isinstance(s, str)})
        return result
    except Exception as e:
        logger.warning(f"Failed to load YAML gene dict: {e}")
        return None


# Built-in fallback dictionary (kept for backward compatibility)
_BUILTIN_GENE_DICT: Dict[str, List[str]] = {
    'coi': list({s.lower() for s in [
        'cytochrome oxidase I',
        'cytochrome oxidase subunit 1 gene',
        'cytochrome c oxidase subunit I',
        'cytochrome oxidase subunit I',
        'cytochrome c oxidase subunit I',
        'cytochrome oxidase subunit 1',
        'mitochondrial COI',
        'cytochrome c oxidase I',
        'cytochrome c oxidase subunit 1',
        'cytochrome c oxidase gene',
        'cox1',
        'cytochrome oxidase subunit I-like',
        'COI', 'CO1', 'COXI', 'COX1'
    ]}),
    '16s': list({s.lower() for s in [
        '16S ribosomal RNA',
        '16S ribosomal RNA gene, partial sequence; mitochondrial',
        'large subunit ribosomal RNA gene, partial sequence; mitochondrial',
        'mitochondrial gene for 16S rRNA, partial sequence',
        '16S ribosomal RNA, partial sequence',
        'large subunit ribosomal RNA',
        '16S rRNA gene', '16S ribosomal RNA gene', '16S-rRNA',
        '16S ribosomal', '16S large subunit rRNA gene',
        'ribosomal RNA large subunit', '16s rRNA', '16S subunit RNA'
    ]}),
    '12s': list({s.lower() for s in [
        '12S ribosomal RNA gene, partial sequence; mitochondrial',
        'mitochondrial gene for 12S ribosomal RNA, patial sequence',
        'small subunit ribosomal RNA gene, partial sequence; mitochondrial',
        'small subunit ribosomal RNA gene',
        'small subunit ribosomal RNA gene, partial sequence; mitochondrial',
        'small subunit ribosomal',
        '12S rRNA gene', '12S rRNA',
        '12S small subunit rRNA gene',
        '12S rRNA, partial sequence',
        'ribosomal RNA small subunit', '12S ribosomal RNA'
    ]}),
    'cob': list({s.lower() for s in [
        'cytochrome b', 'cytb',
        'cytochrome b gene, partial cds; mitochondrial',
        'mitochondrial cytb gene for cytochrome b, partial cds',
        'cytochrome b (Cytb) gene, partial cds; mitochondrial',
        'cytochrome b (CytB) gene, partial cds; mitochondrial',
        'mitochondrial cytb gene for cytochrome b, partial cds',
        'cytochrome oxidase subunit b',
        'COB', 'cytochrome b-like', 'cytB'
    ]}),
    'mtgenome': list({s.lower() for s in [
        'organelle: mitochondrion',
        'mitochondrion, complete genome',
        'mitochondrion, partial genome',
        'mitochondrion sequence',
        'mitochondrial genome',
        'mitochondrial DNA, nearly complete genome',
        'chromosome: MIT',
        'mitochondrial_genome',
        'mt genome', 'mtDNA',
        'complete sequence',
        'complete genome',
        'mitogenome',
        'mitochondrion, complete sequence'
    ]}),
    'cox3': list({s.lower() for s in [
        'cytochrome oxidase subunit III',
        'cytochrome c oxidase subunit 3',
        'cytochrome oxidase 3',
        'cytochrome c oxidase subunit III',
        'cytochrome c oxidase III',
        'cytochrome oxidase subunit 3 mRNA',
        'cytochrome oxidase subunit 3',
        'cytochrome c oxidase subunit III (cox3) gene',
        'COX3', 'CO3', 'COXIII'
    ]}),
    'cox2': list({s.lower() for s in [
        'cytochrome oxidase subunit II',
        'cytochrome oxidase subunit 2',
        'cytochrome c oxidase subunit II',
        'cytochrome c oxidase subunit 2',
        'COX2', 'CO2', 'COXII'
    ]}),
    '18-28s': list({s.lower() for s in [
        '18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S rRNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene',
        '18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S rRNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence',
        'small subunit ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S rRNA gene, and internal transcribed spacer 2, complete sequence; and large subunit ribosomal RNA gene',
        '18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S rRNA gene, internal transcribed spacer 2, and 28S ribosomal RNA gene, complete sequence',
        '18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S rRNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene',
        '18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S rRNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence',
        '18S ribosomal RNA gene, complete sequence; and 28S ribosomal RNA gene, partial sequence'
    ]}),
    'its1-its2': list({s.lower() for s in [
        '5.8S rRNA', '5.8S ribosomal RNA',
        'internal transcribed spacer',  # P0-2a: fixed "trancribed" -> "transcribed"
        'ITS1', 'ITS2',
        'internal transcribed spacer 1 intergenic spacer',
        'internal transcribed spacer 1 intergenic spacer, partial sequence; 5.8S ribosomal RNA gene, complete sequence; and internal transcribed spacer 2 intergenic spacer, partial sequence',
        'ITS'
    ]}),
    '28s': list({s.lower() for s in [
        '28S ribosomal RNA gene',
        '28S ribosomal RNA gene, partial sequence',
        'large subunit ribosomal RNA gene, partial sequence',
        '28S rRNA, partial sequence', '28S rRNA',
        'large subunit ribosomal RNA',
        'ribosomal RNA large subunit gene',
        '28S ribosomal RNA'
    ]}),
    'ef-1': list({s.lower() for s in [
        'elongation factor-1 alpha (EF-1 alpha) mRNA',
        'elongation factor-1', 'elongation factor-1alpha',
        'elongation factor 1-alpha',
        'ef1a', 'ef-1a', 'elongation factor 1 alpha', 'EF1alpha'  # P1-12: added missing synonyms
    ]}),
    '18s': list({s.lower() for s in [
        'small subunit ribosomal RNA gene, partial sequence',
        '18S ribosomal RNA, partial sequence',
        '18s rRNA, complete sequence',
        '18S rRNA, partial sequence',
        '18S ribosomal RNA gene, partial sequence',
        '18S ribosomal RNA, partial sequence',
        '18S ribosomal RNA gene',
        'small subunit ribosomal RNA gene, partial sequence',
        '18S ribosomal RNA', 'small subunit RNA',
        'nuclear small subunit ribosomal RNA',
        'small subunit ribosomal RNA',
        'macronuclear small-subunit ribosomal RNA',
        '18S small subunit ribosomal RNA',
        '18S rRNA', 'small subunit ribosomal RNA gen',
        '18S ribosomal RNA',  # P0-2b: fixed "ribosoaml" -> "ribosomal"
        'ribosomal RNA small-subunit'
    ]}),
    'h3': list({s.lower() for s in [
        'histone H3', 'histone 3', 'histone 3 (H3)',
        'histone 3 (H3) gene, partial cds',
        'histone H3 gene', 'H3 histone',
        'histone (H3) gene', 'histone 3 gene'
    ]})
}

# Initialize from YAML (preferred) or built-in
GENE_DICT: Dict[str, List[str]] = _load_gene_dict_yaml() or _BUILTIN_GENE_DICT

# P2-21: Gene sets are now managed in MatchConfig, not as module-level mutable globals
# We keep these as defaults for init_gene_groups backward compatibility
_DEFAULT_MITO_GENES = {"coi", "cox3", "cox2", "16s", "12s", "cob", "mtgenome"}
_DEFAULT_NUCLEAR_GENES = {"18-28s", "its1-its2", "ef-1", "28s", "18s", "h3"}


# =========================
# P1-9: Pre-compiled regex patterns
# =========================
_COMPILED_PATTERNS: Dict[str, re.Pattern] = {}


def _precompile_patterns():
    """Pre-compile all keyword patterns at module load time."""
    for gene, keywords in GENE_DICT.items():
        for phrase in keywords:
            key = phrase.lower()
            if key in _COMPILED_PATTERNS:
                continue
            escaped = re.escape(key)
            left = r'\b' if re.match(r'\w', key) else r'(?<=\s|^)'
            right = r'\b' if (key and re.match(r'\w', key[-1])) else r'(?=\s|$)'
            _COMPILED_PATTERNS[key] = re.compile(left + escaped + right, re.IGNORECASE)


_precompile_patterns()


def init_gene_groups(mtgenes_list: List[str], ntgenes_list: List[str]) -> Tuple[Set[str], Set[str]]:
    """Initialize mitochondrial/nuclear gene sets. Returns (mito_genes, nuclear_genes)."""
    mito_genes = {g.strip().lower() for g in mtgenes_list}
    nuclear_genes = {g.strip().lower() for g in ntgenes_list}

    for gene in mito_genes:
        if gene not in GENE_DICT:
            logger.warning(f"Mitochondrial gene '{gene}' not in dictionary, adding empty entry")
            GENE_DICT[gene] = []
    for gene in nuclear_genes:
        if gene not in GENE_DICT:
            logger.warning(f"Nuclear gene '{gene}' not in dictionary, adding empty entry")
            GENE_DICT[gene] = []

    # Re-compile after adding new entries
    _precompile_patterns()

    logger.info(f"Mitochondrial genes: {', '.join(sorted(mito_genes))}")
    logger.info(f"Nuclear genes: {', '.join(sorted(nuclear_genes))}")
    return mito_genes, nuclear_genes


# =========================
# Configuration dataclass (P2-21: gene sets now live here)
# =========================
@dataclass
class MatchConfig:
    if_remove_matched: str = "no"

    length2_mtgenes: str = "200:50000"
    length2_ntgenes: str = "200:50000"
    global_length_range: str = "200:50000"

    gene_length_ranges: Dict[str, str] = field(default_factory=lambda: {
        "mtgenome": "3000:500000",
        "18-28s": "5000:15000",
        "coi": "none", "its1-its2": "none", "cob": "none",
        "12s": "none", "16s": "none", "18s": "none",
        "28s": "none", "h3": "none",
    })

    add_assignment_reasons: bool = True
    record_conflict_reasons: bool = True

    # P2-21: Gene sets as config fields instead of mutable globals
    mito_genes: Set[str] = field(default_factory=lambda: _DEFAULT_MITO_GENES.copy())
    nuclear_genes: Set[str] = field(default_factory=lambda: _DEFAULT_NUCLEAR_GENES.copy())

    def get_gene_length_range(self, gene: str) -> str:
        return self.gene_length_ranges.get(gene.lower(), "none")


# =========================
# Load external gene dictionary
# =========================
def load_path_dict(path_dict: str) -> None:
    global GENE_DICT
    if not os.path.exists(path_dict):
        logger.error(f"Gene dictionary file not found: {path_dict}")
        return

    ext = os.path.splitext(path_dict)[1].lower()
    new_dict: Dict[str, List[str]] = {}

    try:
        if ext == ".py":
            spec = importlib.util.spec_from_file_location("ext_dict", path_dict)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            if hasattr(mod, "gene_dict"):
                new_dict = mod.gene_dict
            else:
                logger.error("Python file missing 'gene_dict' variable")
        elif ext == ".txt":
            with open(path_dict, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line and ":" in line:
                        gene, syns = line.split(":", 1)
                        new_dict[gene.strip().lower()] = [
                            s.strip().lower() for s in syns.split(",") if s.strip()
                        ]
        elif ext == ".csv":
            df = pd.read_csv(path_dict, dtype=str)
            if {"gene", "synonyms"}.issubset(df.columns):
                for row in df.itertuples(index=False):
                    gene = row.gene.strip().lower()
                    new_dict[gene] = [
                        s.strip().lower() for s in str(row.synonyms).split(",") if s.strip()
                    ]
            else:
                logger.error("CSV must contain 'gene' and 'synonyms' columns")
        else:
            logger.error(f"Unsupported dictionary format: {ext}")
    except Exception as e:
        logger.error(f"Failed to load external gene dictionary: {e}")

    if new_dict:
        for gene, synonyms in new_dict.items():
            if gene in GENE_DICT:
                GENE_DICT[gene] = list(set(GENE_DICT[gene] + synonyms))
            else:
                GENE_DICT[gene] = synonyms
        _precompile_patterns()
        logger.info(f"External gene dictionary loaded ({len(new_dict)} genes)")



def smart_match(phrase: str, text: str) -> bool:
    """Smart phrase matching using pre-compiled patterns (P1-9)."""
    key = phrase.lower()
    pat = _COMPILED_PATTERNS.get(key)
    if pat is None:
        # Fallback for phrases not in GENE_DICT (e.g. from external dicts)
        escaped = re.escape(key)
        left = r'\b' if re.match(r'\w', key) else r'(?<=\s|^)'
        right = r'\b' if (key and re.match(r'\w', key[-1])) else r'(?=\s|$)'
        pat = re.compile(left + escaped + right, re.IGNORECASE)
        _COMPILED_PATTERNS[key] = pat
    return pat.search(text) is not None


def any_match(word_list: List[str], text: str) -> bool:
    return any(smart_match(w, text) for w in word_list)


def clean_length(length_str: str) -> int:
    try:
        return int(re.sub(r'[^\d]', '', str(length_str)))
    except (ValueError, TypeError):
        return 0


def filter_dataframe(df: pd.DataFrame, moleculetype: str, mol_type: str,
                     length_filter: str, organelle_filter: str) -> pd.DataFrame:
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


DIAG_COLUMNS = ["gene_type", "Conflict", "match_source", "Original_match",
                "Conflict_reason", "Assignment_reason"]


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    front = [c for c in DIAG_COLUMNS if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest]


# =========================
# P1-8: Shared gene matching loop
# =========================
def _match_gene_loop(
    definition: str,
    length_val: int,
    config: MatchConfig,
    gene_set: Set[str],
    skip_genes: Set[str],
    priority_genes: List[str],
    global_range: str,
    extra_priority_check=None,
) -> Tuple[Set[str], List[str]]:
    """Shared matching loop for both mitochondrial and nuclear gene sets."""
    matched = set()
    notes = []

    if not length_in_range(length_val, global_range):
        notes.append(f"Length {length_val} outside range {global_range}")
        return matched, notes

    # Priority genes
    for pg in priority_genes:
        if pg in skip_genes or pg not in gene_set:
            continue
        if any_match(GENE_DICT.get(pg, []), definition):
            range_str = config.get_gene_length_range(pg)
            if length_in_range(length_val, range_str):
                # Extra check for mtgenome
                if pg == "mtgenome" and length_val <= 3000:
                    notes.append(f"mtgenome: keyword matched but Length={length_val}<=3000")
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
        keywords = GENE_DICT.get(gene, [])
        if not keywords:
            continue
        if not any_match(keywords, definition):
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


def match_mito_genes(definition: str, length_val: int, molecule_type: str,
                     config: MatchConfig) -> Tuple[Set[str], List[str]]:
    """Mitochondrial gene matching: mtgenome priority, then other mito genes."""
    return _match_gene_loop(
        definition, length_val, config,
        gene_set=config.mito_genes,
        skip_genes=set(),
        priority_genes=["mtgenome"],
        global_range=config.length2_mtgenes,
    )


def match_nuclear_genes(definition: str, length_val: int,
                        config: MatchConfig) -> Tuple[Set[str], List[str]]:
    """Nuclear gene matching: 18-28s priority, then its1-its2, then other nuclear genes."""
    return _match_gene_loop(
        definition, length_val, config,
        gene_set=config.nuclear_genes,
        skip_genes=set(),
        priority_genes=["18-28s", "its1-its2"],
        global_range=config.length2_ntgenes,
    )


def first_round_match(row: Dict[str, Any], config: MatchConfig) -> Dict[str, Any]:
    definition = str(row.get("Definition", "")).lower().strip()
    organelle = str(row.get("organelle", "")).lower().strip()
    length_val = clean_length(row.get("Length", "0"))

    group = "mito" if "mitochondrion" in organelle else "nuclear"
    logger.info(f"Round 1 LocusID={row.get('LocusID', 'NA')}, group={group}")
    logger.debug(f"  Definition: {definition[:120]}")

    if group == "mito":
        molecule_type = str(row.get("MoleculeType", "")).lower().strip()
        matched, notes = match_mito_genes(definition, length_val, molecule_type, config)
    else:
        # P2-19: Removed unreachable dead code (mito check in nuclear branch)
        matched, notes = match_nuclear_genes(definition, length_val, config)

    gene_type = ",".join(sorted(matched))
    match_source = "first_match" if matched else ""
    return _build_result(gene_type, match_source, notes)


# =========================
# Round 1.5: Conflict check
# =========================
def _check_ribosomal_conflict(gene: str, definition: str, length_val: int) -> Optional[Dict]:
    if gene == "12s" and re.search(r'large subunit ribosomal', definition, re.IGNORECASE):
        return {
            "resolved": {"12s"},  # P0-4: preserve original (consistent with 18s/28s now)
            "reason": "12s matched with 'large subunit ribosomal' in definition, possible 16S mismatch, flagged",
            "conflict_reason": "12s conflicts with 'large subunit ribosomal' description"
        }
    if gene == "16s" and re.search(r'small subunit ribosomal', definition, re.IGNORECASE):
        return {
            "resolved": {"16s"},
            "reason": "16s matched with 'small subunit ribosomal' in definition, possible 12S mismatch, flagged",
            "conflict_reason": "16s conflicts with 'small subunit ribosomal' description"
        }
    if gene in ("18s", "28s") and "mitochondrial" in definition.lower():
        # P0-4: Always preserve original gene name, record conflict in metadata
        if length_val > 800:
            return {
                "resolved": {gene},
                "reason": f"{gene} in definition with 'mitochondrial' but length {length_val}>800, kept {gene}",
                "conflict_reason": f"{gene} conflicts with mitochondrial annotation"
            }
        else:
            return {
                "resolved": {gene},  # P0-4: no more '?' suffix — consistent with 12s/16s
                "reason": f"{gene} in definition with 'mitochondrial' and length {length_val}<=800, flagged as uncertain",
                "conflict_reason": f"{gene} conflicts with mitochondrial annotation, short length"
            }
    return None


def conflict_check_and_reassign(gene_type_str: str, definition: str,
                                topology: str, group: str,
                                length_val: int) -> Dict[str, Any]:
    genes = {g.strip().lower() for g in gene_type_str.split(",") if g.strip()}
    original = gene_type_str
    notes = []

    skip_types = {"mtgenome", "18-28s", "its1-its2"}
    if genes & skip_types:
        notes.append("Special type (mtgenome/18-28s/its1-its2), skip conflict check")
        return _build_result(gene_type_str, "", notes, original=original)

    conflict = False
    conflict_reason = ""
    resolved = set(genes)

    if group.lower() == "mito":
        resolved, conflict, conflict_reason, notes = _conflict_mito(
            genes, definition, length_val, notes)
    else:
        resolved, conflict, conflict_reason, notes = _conflict_nuclear(
            genes, definition, length_val, notes)

    new_gene_type = ",".join(sorted(resolved))

    if len(resolved) > 1 and not conflict:
        conflict = True
        conflict_reason = f"Multiple candidate gene types: {', '.join(sorted(resolved))}"
        notes.append("Multiple candidates flagged as conflict")

    match_source = "reassign_mode" if new_gene_type != original else ""
    if new_gene_type != original:
        notes.append(f"Reassigned: [{original}] -> [{new_gene_type}]")

    return _build_result(new_gene_type, match_source, notes,
                         original=original, conflict=conflict,
                         conflict_reason=conflict_reason)


def _conflict_mito(genes: Set[str], definition: str, length_val: int,
                   notes: List[str]) -> Tuple[Set[str], bool, str, List[str]]:
    conflict = False
    conflict_reason = ""
    resolved = set(genes)

    if len(genes) > 2:
        if length_val > 3000:
            resolved = {"mtgenome"}
            notes.append(f"Mitochondrial multi-match ({len(genes)}) + Length>3000 -> reassigned to mtgenome")
        else:
            conflict = True
            conflict_reason = f"Mitochondrial multi-match ({len(genes)}) but length {length_val}<3000"
            notes.append("Mitochondrial multi-match but insufficient length, kept original + flagged")
        return resolved, conflict, conflict_reason, notes

    if len(genes) == 2:
        has_trna = re.search(r'trna', definition, re.IGNORECASE)
        has_complete = re.search(r'complete sequence', definition, re.IGNORECASE)
        if has_trna or has_complete:
            conflict = True
            trigger = "tRNA" if has_trna else "complete sequence"
            conflict_reason = f"Two mitochondrial genes + definition contains '{trigger}'"
            notes.append(f"Dual match + '{trigger}' in definition, kept original + flagged")
            return resolved, conflict, conflict_reason, notes

    for gene in list(genes):
        result = _check_ribosomal_conflict(gene, definition, length_val)
        if result:
            resolved = result["resolved"]
            conflict = True
            conflict_reason = result["conflict_reason"]
            notes.append(result["reason"])
            return resolved, conflict, conflict_reason, notes

    return resolved, conflict, conflict_reason, notes


def _conflict_nuclear(genes: Set[str], definition: str, length_val: int,
                      notes: List[str]) -> Tuple[Set[str], bool, str, List[str]]:
    conflict = False
    conflict_reason = ""
    resolved = set(genes)

    for gene in list(genes):
        result = _check_ribosomal_conflict(gene, definition, length_val)
        if result:
            resolved = result["resolved"]
            conflict = True
            conflict_reason = result["conflict_reason"]
            notes.append(result["reason"])
            return resolved, conflict, conflict_reason, notes

    return resolved, conflict, conflict_reason, notes


# =========================
# Round 2: Recheck unmatched
# =========================
def recheck_match_row(row: Dict[str, Any], config: MatchConfig) -> Dict[str, Any]:
    definition = str(row.get("Definition", "")).lower().strip()
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

    if topology == "circular":
        candidates.add("mtgenome")
        notes.append("Topology=circular -> direct mtgenome")
    else:
        for gene, keywords in GENE_DICT.items():
            if not keywords:
                continue
            if not any_match(keywords, definition):
                continue

            if gene == "mtgenome":
                if (length_val > 3000 and molecule_type == "dna"
                        and (has_mito_def or has_mito_org)):
                    candidates.add("mtgenome")
                    notes.append("mtgenome: keyword + Length>3000 + DNA + mitochondrial features")
                else:
                    conflict_notes.append(
                        f"mtgenome: conditions not met (Length={length_val}, "
                        f"MoleculeType={molecule_type}, "
                        f"mito_in_def={has_mito_def}, mito_in_org={has_mito_org})")
            else:
                candidates.add(gene)
                notes.append(f"{gene}: keyword matched")

        if not candidates and pcr_primers:
            coi_primers = ["HCO2198", "lc01490", "LCO1490"]
            if any(re.search(rf'(?:fwd|rev)_name:\s*{re.escape(p)}', pcr_primers, re.IGNORECASE)
                   for p in coi_primers):
                candidates.add("coi")
                notes.append("PCR primer match (HCO2198/lc01490) -> coi")

            if re.search(r'(?:fwd|rev)_name:\s*28SC1', pcr_primers, re.IGNORECASE):
                candidates.add("28s")
                notes.append("PCR primer match (28SC1) -> 28s")

    if candidates:
        if len(candidates) > 1:
            conflict_notes.append(f"Multiple candidates: {', '.join(sorted(candidates))}")
        elif len(candidates) == 1:
            gene_assigned = next(iter(candidates))
            if gene_assigned in config.mito_genes and not (has_mito_def or has_mito_org):
                conflict_notes.append(
                    f"Mitochondrial candidate ({gene_assigned}) but no mitochondrial marker in organelle/definition")
            if gene_assigned in config.nuclear_genes and (has_mito_def or has_mito_org):
                conflict_notes.append(
                    f"Nuclear candidate ({gene_assigned}) but mitochondrial marker in organelle/definition")

    gene_type = ",".join(sorted(candidates))
    match_source = "recheck_mode" if gene_type else ""
    conflict = bool(conflict_notes) and bool(gene_type)

    return _build_result(
        gene_type, match_source, notes,
        original=gene_type,
        conflict=conflict,
        conflict_reason="; ".join(conflict_notes) if conflict else ""
    )


# =========================
# Result builder
# =========================
def _build_result(gene_type: str, match_source: str, notes: List[str],
                  original: str = "", conflict: bool = False,
                  conflict_reason: str = "") -> Dict[str, Any]:
    return {
        "gene_type": gene_type,
        "match_source": match_source,
        "Original_match": original or gene_type,
        "conflict": conflict,
        "Conflict_reason": conflict_reason,
        "Assignment_reason": "; ".join(notes)
    }


# =========================
# Row processing
# =========================
def process_row(row: Dict[str, Any], config: MatchConfig) -> Dict[str, Any]:
    prelim = first_round_match(row, config)
    gene_type_str = prelim["gene_type"]

    if not gene_type_str:
        return prelim

    definition = str(row.get("Definition", "")).lower()
    group = "mito" if "mitochondrion" in str(row.get("organelle", "")).lower() else "nuclear"
    length_val = clean_length(row.get("Length", "0"))

    final = conflict_check_and_reassign(gene_type_str, definition, "", group, length_val)

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
        "Assignment_reason": "; ".join(reasons)
    }


# =========================
# Main flows
# =========================
def process_prematch(input_file: str, output_dir: str, config: MatchConfig,
                     prematch_output: str, multip_output: str,
                     unmatched_output: str, conflicted_output: str,
                     moleculetype: str, mol_type: str,
                     length_filter: str, organelle_filter: str) -> pd.DataFrame:
    os.makedirs(output_dir, exist_ok=True)

    ext = os.path.splitext(input_file)[1].lower()
    read_funcs = {
        ".csv": lambda f: pd.read_csv(f, dtype=str),
        ".tsv": lambda f: pd.read_csv(f, sep="\t", dtype=str),
        ".xls": lambda f: pd.read_excel(f, dtype=str),
        ".xlsx": lambda f: pd.read_excel(f, dtype=str),
    }
    if ext not in read_funcs:
        logger.error(f"Unsupported format: {ext}")
        sys.exit(1)
    df = read_funcs[ext](input_file)
    df.columns = df.columns.str.strip()

    df["Length"] = df["Length"].astype(str).str.replace(r'\s*bp\s*', '', regex=True)
    df["Length"] = pd.to_numeric(df["Length"], errors='coerce')

    lo, hi = parse_interval(config.global_length_range)
    if lo is not None:
        df = df[df["Length"] >= lo]
    if hi is not None:
        df = df[df["Length"] <= hi]
    logger.info(f"After global length filter: {len(df)} records")

    df = filter_dataframe(df, moleculetype, mol_type, length_filter, organelle_filter)
    logger.info(f"After additional filters: {len(df)} records")

    results = []
    conflict_records = []
    for row in df.itertuples(index=False):
        row_dict = row._asdict()
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

    res_df = pd.DataFrame(results)
    if "LocusID" in res_df.columns:
        res_df = res_df.drop_duplicates(subset=["LocusID"])

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
        conflict_df = reorder_columns(pd.DataFrame(conflict_records))
        _save_csv(conflict_df, output_dir, conflicted_output, "Conflict")

    multip_df = assigned_df[assigned_df["gene_type"].apply(
        lambda x: len(x.split(",")) > 1 if x else False)]
    if not multip_df.empty:
        _save_csv(multip_df, output_dir, multip_output, "Multi-match")

    total = len(res_df)
    logger.info(f"Round 1 summary: total={total}, matched={len(assigned_df)}, "
                f"conflict={len(conflict_records)}, unmatched={len(unmatched_df)}")
    return res_df


def process_recheck(input_file: str, output_dir: str, config: MatchConfig,
                    apply_length_filter: bool = True,
                    assigned_out: str = "assigned_genes_types2.csv",
                    multip_out: str = "multiple_matched2.csv",
                    unmatched_out: str = "unmatched_sequences2.csv") -> pd.DataFrame:
    os.makedirs(output_dir, exist_ok=True)

    ext = os.path.splitext(input_file)[1].lower()
    read_funcs = {
        ".csv": lambda f: pd.read_csv(f, dtype=str),
        ".tsv": lambda f: pd.read_csv(f, sep="\t", dtype=str),
        ".xls": lambda f: pd.read_excel(f, dtype=str),
        ".xlsx": lambda f: pd.read_excel(f, dtype=str),
    }
    if ext not in read_funcs:
        logger.error(f"Unsupported format: {ext}")
        sys.exit(1)
    df = read_funcs[ext](input_file)
    df.columns = df.columns.str.strip()

    if df.empty:
        logger.info("Round 2: no unmatched records to recheck")
        _save_csv(df, output_dir, assigned_out, "Round 2 matched")
        _save_csv(df, output_dir, unmatched_out, "Round 2 unmatched")
        return df

    # P0-3: Use union of mt+nt ranges instead of nuclear-only
    if apply_length_filter:
        df["Length"] = df["Length"].astype(str).str.replace(r'\s*bp\s*', '', regex=True)
        df["Length"] = pd.to_numeric(df["Length"], errors='coerce')
        lo_mt, hi_mt = parse_interval(config.length2_mtgenes)
        lo_nt, hi_nt = parse_interval(config.length2_ntgenes)
        lo = min(lo_mt or 0, lo_nt or 0)
        hi = max(hi_mt or float('inf'), hi_nt or float('inf'))
        if lo is not None:
            df = df[df["Length"] >= lo]
        if hi is not None and hi != float('inf'):
            df = df[df["Length"] <= hi]

    results = []
    conflict_records = []
    for row in df.itertuples(index=False):
        row_dict = row._asdict()
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

    res_df = pd.DataFrame(results)
    if "LocusID" in res_df.columns:
        res_df = res_df.drop_duplicates(subset=["LocusID"])

    if not config.add_assignment_reasons and "Assignment_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Assignment_reason"])
    if not config.record_conflict_reasons and "Conflict_reason" in res_df.columns:
        res_df = res_df.drop(columns=["Conflict_reason"])

    res_df = reorder_columns(res_df)

    assigned_df = res_df[res_df["gene_type"] != ""]
    unmatched_df = res_df[res_df["gene_type"] == ""]
    multi_df = assigned_df[assigned_df["gene_type"].apply(
        lambda x: len(x.split(",")) > 1 if x else False)]

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

    logger.info(f"Round 2 summary: total={len(res_df)}, matched={len(assigned_df)}, "
                f"conflict={len(conflict_records)}, unmatched={len(unmatched_df)}")
    return res_df


def extract_genes(df: pd.DataFrame, gene_types: Optional[List[str]],
                  output_dir: str) -> None:
    if not gene_types:
        all_genes: Set[str] = set()
        for gt in df["gene_type"].dropna():
            all_genes.update(s.strip() for s in str(gt).split(",") if s.strip())
        gene_types = sorted(all_genes)
        logger.info(f"Auto-detected gene types: {gene_types}")

    for gene in gene_types:
        subset = df[df["gene_type"].apply(
            lambda x: gene in [s.strip() for s in str(x).split(",")] if x else False)]
        if subset.empty:
            logger.warning(f"No records for gene type '{gene}'")
            continue
        safe_name = re.sub(r'[^\w\-]', '_', gene)
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
    if 'whole genome shotgun' in def_lower:
        return 'whole_genome_shotgun'
    elif 'genome assembly' in def_lower:
        return 'genome_assembly'
    elif 'tpa_asm' in def_lower or 'tpa:' in def_lower:
        return 'TPA_assembly'
    elif 'chromosome:' in def_lower:
        return 'chromosome_data'
    elif 'scaffold' in def_lower:
        return 'scaffold'
    elif 'contig' in def_lower:
        return 'contig'

    # RNA types
    elif 'mrna sequence' in def_lower or 'cdna clone' in def_lower:
        return 'mRNA/cDNA'
    elif 'trna' in def_lower:
        return 'tRNA'
    elif 'rrna' in def_lower and 'mitochondrial' not in def_lower:
        return 'nuclear_rRNA'

    # Other genes not in standard list
    elif 'nad' in def_lower or 'nadh' in def_lower:
        return 'NADH_dehydrogenase'
    elif 'atp synthase' in def_lower or 'atpase' in def_lower:
        return 'ATP_synthase'
    elif 'creatine kinase' in def_lower:
        return 'creatine_kinase'
    elif 'cytochrome b' in def_lower or 'cytb' in def_lower:
        return 'cytochrome_b_variant'
    elif 'nadh' in def_lower:
        return 'NADH_dehydrogenase'
    else:
        return 'other'


def generate_unmatched_report(unmatched_df: pd.DataFrame, output_dir: str) -> str:
    """Generate a summary report of unmatched sequences by category."""
    if unmatched_df.empty:
        return ""

    # Categorize each unmatched sequence
    unmatched_df = unmatched_df.copy()
    unmatched_df['category'] = unmatched_df['Definition'].apply(categorize_unmatched_definition)

    # Generate summary
    summary_lines = []
    summary_lines.append("=" * 60)
    summary_lines.append("UNMATCHED SEQUENCES REPORT")
    summary_lines.append("=" * 60)
    summary_lines.append(f"Total unmatched: {len(unmatched_df)}")
    summary_lines.append("")

    category_counts = unmatched_df['category'].value_counts()
    summary_lines.append("Category breakdown:")
    summary_lines.append("-" * 40)

    for cat, count in category_counts.items():
        pct = count / len(unmatched_df) * 100
        summary_lines.append(f"  {cat}: {count} ({pct:.1f}%)")

    summary_lines.append("")
    summary_lines.append("Sample definitions per category:")
    summary_lines.append("-" * 40)

    for cat in category_counts.index[:5]:  # Top 5 categories
        samples = unmatched_df[unmatched_df['category'] == cat]['Definition'].head(3).tolist()
        summary_lines.append(f"\n[{cat}]")
        for s in samples:
            summary_lines.append(f"  - {str(s)[:100]}...")

    summary_lines.append("")
    summary_lines.append("=" * 60)

    # Save report
    report_path = os.path.join(output_dir, "unmatched_report.txt")
    report_content = "\n".join(summary_lines)
    with open(report_path, 'w', encoding='utf-8') as f:
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
    logger.info("=" * 60)
    logger.info("Gene Type Classifier (refactored)")
    logger.info("=" * 60)

    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    mtgenes = [s.strip() for s in args.mtgenes_list.split(",") if s.strip()]
    ntgenes = [s.strip() for s in args.ntgenes_list.split(",") if s.strip()]
    mito_genes, nuclear_genes = init_gene_groups(mtgenes, ntgenes)

    if args.path_dict:
        load_path_dict(args.path_dict)

    gene_ranges = {
        "mtgenome": args.length_range2_mtgenome,
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
        global_length_range=args.length_range2_all,
        gene_length_ranges=gene_ranges,
        add_assignment_reasons=args.if_add_assigned_gene_types_reasons.lower() == "yes",
        record_conflict_reasons=args.if_recorded_conflict_reasons.lower() == "yes",
        mito_genes=mito_genes,
        nuclear_genes=nuclear_genes,
    )

    # Round 1
    logger.info("-" * 40)
    logger.info("Round 1: Gene type matching")
    logger.info("-" * 40)
    prematch_df = process_prematch(
        args.input, output_dir, config,
        args.prematch_output_file_name,
        args.multip_matched_genestype,
        args.unmatched_sequences,
        args.conflicted_sequences,
        args.moleculetype, args.mol_type,
        args.length, args.organelle
    )

    # Gene extraction
    gene_types_arg = args.which_gene_types_extract
    if gene_types_arg:
        extract_genes(prematch_df,
                      [s.strip().lower() for s in gene_types_arg.split(",")],
                      output_dir)

    # Round 2
    if args.if_recheck_unmatched.lower() == "yes":
        unmatched_file = os.path.join(output_dir, args.unmatched_sequences)
        if os.path.exists(unmatched_file) and os.path.getsize(unmatched_file) > 50:
            logger.info("-" * 40)
            logger.info("Round 2: Unmatched record recheck")
            logger.info("-" * 40)

            process_recheck(
                unmatched_file, output_dir, config,
                apply_length_filter=(args.if_length_range_work2_recheck.lower() == "yes")
            )

            round1_output_path = os.path.join(output_dir, args.prematch_output_file_name)
            round2_output_path = os.path.join(output_dir, "assigned_genes_types2.csv")
            if os.path.exists(round1_output_path) and os.path.exists(round2_output_path):
                try:
                    round1_df = pd.read_csv(round1_output_path, dtype=str)
                    round2_df = pd.read_csv(round2_output_path, dtype=str)
                    merged = pd.concat([round1_df, round2_df])
                    if "LocusID" in merged.columns:
                        merged = merged.drop_duplicates(subset=["LocusID"])
                    merged = reorder_columns(merged)
                    all_file = os.path.join(output_dir, "assigned_genes_types_all.csv")
                    merged.to_csv(all_file, index=False)
                    logger.info(f"Merged results: {all_file} ({len(merged)} records)")
                except Exception as e:
                    logger.error(f"Merge failed: {e}")
        else:
            logger.warning(f"Unmatched file not found: {unmatched_file}, skipping round 2")

    logger.info("=" * 60)
    logger.info("Processing complete")
    logger.info("=" * 60)


def classify(
    input_file: str,
    output_dir: str,
    config: MatchConfig = None,
    prematch_output: str = "assigned_genes_type.csv",
    unmatched_output: str = "unmatched_sequences.csv",
    if_recheck: bool = True,
    apply_length_filter_recheck: bool = True,
    extract_gene_types: list = None,
) -> StepResult:
    """Classify gene types from metadata CSV. Returns StepResult."""
    import time
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    if config is None:
        config = MatchConfig()

    try:
        prematch_df = process_prematch(
            input_file, output_dir, config,
            prematch_output, "multiple_matched.csv",
            unmatched_output, "conflicted_genes.csv",
            "all", "all", "all", "all",
        )

        if if_recheck:
            unmatched_file = os.path.join(output_dir, unmatched_output)
            if os.path.exists(unmatched_file) and os.path.getsize(unmatched_file) > 50:
                process_recheck(
                    unmatched_file, output_dir, config,
                    apply_length_filter=apply_length_filter_recheck,
                )
                round1_output_path = os.path.join(output_dir, prematch_output)
                round2_output_path = os.path.join(output_dir, "assigned_genes_types2.csv")
                if os.path.exists(round1_output_path) and os.path.exists(round2_output_path):
                    try:
                        round1_df = pd.read_csv(round1_output_path, dtype=str)
                        round2_df = pd.read_csv(round2_output_path, dtype=str)
                        merged = pd.concat([round1_df, round2_df])
                        if "LocusID" in merged.columns:
                            merged = merged.drop_duplicates(subset=["LocusID"])
                        merged = reorder_columns(merged)
                        all_file = os.path.join(output_dir, "assigned_genes_types_all.csv")
                        merged.to_csv(all_file, index=False)
                    except Exception:
                        pass

        if extract_gene_types:
            extract_genes(prematch_df, extract_gene_types, output_dir)

        # Generate unmatched report
        unmatched_file = os.path.join(output_dir, unmatched_output)
        if os.path.exists(unmatched_file):
            try:
                unmatched_df = pd.read_csv(unmatched_file, dtype=str)
                if not unmatched_df.empty:
                    generate_unmatched_report(unmatched_df, output_dir)
            except Exception as e:
                logger.warning(f"Failed to generate unmatched report: {e}")

        all_file = os.path.join(output_dir, "assigned_genes_types_all.csv")
        if os.path.exists(all_file):
            rows = len(pd.read_csv(all_file, dtype=str))
        else:
            rows = 0

        return StepResult(
            success=True,
            output_file=all_file,
            rows=rows,
            elapsed=time.time() - start_time,
        )
    except Exception as e:
        logger.error(f"Classify failed: {e}")
        return StepResult(success=False, output_file="", elapsed=time.time() - start_time)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Gene Type Classifier: GenBank metadata -> gene_type labels",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True, help='Input file (CSV/TSV/Excel)')
    parser.add_argument('-o', '--output', default='assign_genes_type', help='Output directory')

    parser.add_argument('--prematch_output_file_name', default='assigned_genes_type.csv')
    parser.add_argument('--multip_matched_genestype', default='multiple_matched.csv')
    parser.add_argument('--unmatched_sequences', default='unmatched_sequences.csv')
    parser.add_argument('--conflicted_sequences', default='conflicted_genes.csv')

    parser.add_argument('--if_remove_matched_sequences', '-ifrms', default='no')

    parser.add_argument('--length_range2_mtgenome', default='3000:500000')
    parser.add_argument('--length_range2_18_28s', default='5000:15000')
    parser.add_argument('--length_range2_coi', default='none')
    parser.add_argument('--length_range2_its', default='none')
    parser.add_argument('--length_range2_cob', default='none')
    parser.add_argument('--length_range2_12s', default='none')
    parser.add_argument('--length_range2_16s', default='none')
    parser.add_argument('--length_range2_18s', default='none')
    parser.add_argument('--length_range2_28s', default='none')
    parser.add_argument('--length_range2_h3', default='none')

    parser.add_argument('--gene_length_range', nargs='*', default=None,
                        help='Per-gene length ranges: gene=lower:upper')

    parser.add_argument('--length2_mtgenes', default='200:50000')
    parser.add_argument('--length2_ntgenes', default='200:50000')
    parser.add_argument('--length_range2_all', default='200:50000')
    parser.add_argument('--length', default='all')
    parser.add_argument('--organelle', default='all')
    parser.add_argument('--moleculetype', default='all')
    parser.add_argument('--mol_type', default='all')

    parser.add_argument('--mtgenes_list', default='coi,cox3,cox2,16s,12s,cob,mtgenome')
    parser.add_argument('--ntgenes_list', default='18-28s,its1-its2,ef-1,28s,18s,h3')
    parser.add_argument('--path_dict', default=None, help='External gene dictionary file')

    # P2-15: Fixed typo, with backward-compatible alias
    parser.add_argument('--which_gene_types_extract', '--which_gene_types_extact', '-wgte', default=None)

    parser.add_argument('--if_recheck_unmatched', default='yes')
    parser.add_argument('--if_length_range_work2_recheck', default='yes')
    parser.add_argument('--if_add_assigned_gene_types_reasons', default='yes')
    parser.add_argument('--if_recorded_conflict_reasons', default='yes')

    args = parser.parse_args(argv)

    check_argv = argv if argv else sys.argv
    if '--which_gene_types_extact' in check_argv:
        warnings.warn("--which_gene_types_extact is deprecated (typo), use --which_gene_types_extract",
                       DeprecationWarning, stacklevel=2)

    _run_cli(args)


if __name__ == "__main__":
    main()
