#!/usr/bin/env python3
"""
GenBank Metadata Extractor (refactored)
========================================
Refactored from gb_baseinformation_extractv16.7.py with:
  - Shared _process_record_loop eliminates triple code duplication (P1-7)
  - Shared taxonomy cache via multiprocessing.Manager (P1-10)
  - Combined validate + count in single file read (P1-11)
  - Renamed parameter to avoid shadowing (P1-20)
  - Removed duplicate format_size/get_file_size_str
  - Removed unused _ACCESSION_RE
"""

import os
import re
import sys
import csv
import json
import shutil
import argparse
import traceback
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Callable
from g2t.utils import StepResult
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import concurrent.futures
import multiprocessing
import importlib.util
import time
import logging
from enum import Enum

# ===============================
# Logging
# ===============================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    encoding='utf-8'
)
logger = logging.getLogger(__name__)

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

# ===============================
# Global state
# ===============================
ete3_available = False
NCBITaxa = None
TAXONOMY_CACHE: Dict[str, Dict[str, str]] = {}
TAXONOMY_ERROR_REPORTED: set = set()

try:
    import psutil
except ImportError:
    psutil = None

GB_EXTENSIONS = (".gb", ".gbk", ".genbank", ".gbff")

# ===============================
# Assembly regex (pre-compiled)
# ===============================
_FLEX_AM = re.compile(r'Assembly\s+Method\s*::?\s*(.*?)(?:\n|$)', re.IGNORECASE)
_FLEX_AN = re.compile(r'Assembly\s+Name\s*::?\s*(.*?)(?:\n|$)', re.IGNORECASE)
_FLEX_ST = re.compile(r'Sequencing\s+Technology\s*::?\s*(.*?)(?:\n|$)', re.IGNORECASE)
_FLEX_GC = re.compile(r'Genome\s+Coverage\s*::?\s*(.*?)(?:\n|$)', re.IGNORECASE)
_FLEX_EL = re.compile(r'Expected\s+Final\s+Version\s*::?\s*(.*?)(?:\n|$)', re.IGNORECASE)
_BLOCK_RE = re.compile(
    r'##Assembly-Data-START##(.*?)##Assembly-Data-END##', re.DOTALL
)

# ===============================
# Data classes
# ===============================
class FileStatus(Enum):
    SUCCESS = "success"
    PARTIAL = "partial"
    FAILED = "failed"
    SKIPPED = "skipped"
    INVALID = "invalid"
    EMPTY = "empty"


@dataclass
class RecordError:
    file_path: str
    record_index: int
    accession: str = ""
    error_type: str = ""
    error_message: str = ""
    phase: str = ""


@dataclass
class RecordTrack:
    accession: str
    record_index: int
    metadata_ok: bool = False
    assembly_ok: bool = False
    assembly_source: str = ""
    has_assembly_data: bool = False
    error_message: str = ""


@dataclass
class FileResult:
    file_path: str
    status: FileStatus = FileStatus.FAILED
    total_records: int = 0
    metadata_extracted: int = 0
    metadata_failed: int = 0
    assembly_extracted: int = 0
    assembly_no_data: int = 0
    assembly_failed: int = 0
    error_message: str = ""
    elapsed_time: float = 0.0
    record_errors: List[RecordError] = field(default_factory=list)
    record_tracks: List[RecordTrack] = field(default_factory=list)

    def add_record_error(self, error: RecordError):
        self.record_errors.append(error)

    @property
    def extracted_records(self) -> int:
        return self.metadata_extracted

    @property
    def failed_records(self) -> int:
        return self.metadata_failed


@dataclass
class ProcessingSummary:
    total_files_found: int = 0
    gb_files_found: int = 0
    non_gb_files_skipped: int = 0
    files_succeeded: int = 0
    files_partial: int = 0
    files_failed: int = 0
    files_empty: int = 0
    files_invalid: int = 0
    total_records_in_files: int = 0
    total_metadata_extracted: int = 0
    total_metadata_failed: int = 0
    total_assembly_extracted: int = 0
    total_assembly_no_data: int = 0
    total_assembly_failed: int = 0
    total_elapsed_time: float = 0.0
    file_results: List[FileResult] = field(default_factory=list)

    def add_result(self, result: FileResult):
        self.file_results.append(result)
        _map = {
            FileStatus.SUCCESS: 'files_succeeded',
            FileStatus.PARTIAL: 'files_partial',
            FileStatus.FAILED: 'files_failed',
            FileStatus.EMPTY: 'files_empty',
            FileStatus.INVALID: 'files_invalid',
        }
        attr = _map.get(result.status)
        if attr:
            setattr(self, attr, getattr(self, attr) + 1)
        self.total_records_in_files += result.total_records
        self.total_metadata_extracted += result.metadata_extracted
        self.total_metadata_failed += result.metadata_failed
        self.total_assembly_extracted += result.assembly_extracted
        self.total_assembly_no_data += result.assembly_no_data
        self.total_assembly_failed += result.assembly_failed

    def get_failed_files(self) -> List[str]:
        return [r.file_path for r in self.file_results
                if r.status in (FileStatus.FAILED, FileStatus.INVALID, FileStatus.EMPTY)]

    def get_all_record_errors(self) -> List[RecordError]:
        errs = []
        for r in self.file_results:
            errs.extend(r.record_errors)
        return errs

    def print_summary(self):
        logger.info("")
        logger.info("=" * 75)
        logger.info("                          Summary")
        logger.info("=" * 75)
        logger.info(f"  Total files:              {self.total_files_found}")
        logger.info(f"    GenBank files:          {self.gb_files_found}")
        logger.info(f"    Non-GB (skipped):       {self.non_gb_files_skipped}")
        logger.info("-" * 75)
        logger.info("  File results:")
        logger.info(f"    Success:               {self.files_succeeded}")
        logger.info(f"    Partial:               {self.files_partial}")
        logger.info(f"    Failed:                {self.files_failed}")
        logger.info(f"    Empty:                 {self.files_empty}")
        logger.info(f"    Invalid:               {self.files_invalid}")
        logger.info("-" * 75)
        logger.info("  Record statistics:")
        logger.info(f"    Total records:         {self.total_records_in_files}")
        logger.info(f"    Metadata extracted:    {self.total_metadata_extracted}")
        logger.info(f"    Metadata failed:       {self.total_metadata_failed}")
        md_rate = (self.total_metadata_extracted / self.total_records_in_files * 100
                   if self.total_records_in_files > 0 else 0)
        logger.info(f"    Metadata rate:         {md_rate:.1f}%")
        logger.info("")
        logger.info(f"    Assembly extracted:    {self.total_assembly_extracted}")
        logger.info(f"    Assembly no data:      {self.total_assembly_no_data}")
        logger.info(f"    Assembly failed:       {self.total_assembly_failed}")
        asm_with_data = self.total_assembly_extracted + self.total_assembly_failed
        asm_rate = (self.total_assembly_extracted / asm_with_data * 100
                    if asm_with_data > 0 else 0)
        logger.info(f"    Assembly coverage:     {self.total_assembly_extracted}/{asm_with_data} ({asm_rate:.1f}%)")
        logger.info("-" * 75)
        logger.info(f"  Total time:              {self.total_elapsed_time:.2f}s")
        logger.info("=" * 75)
        self._print_per_file_details()

    def _print_per_file_details(self):
        if not self.file_results:
            return
        logger.info("")
        logger.info("  Per-file details:")
        logger.info("  " + "-" * 73)
        for r in self.file_results:
            fname = os.path.basename(r.file_path)
            icon = {"success": "OK", "partial": "~", "failed": "X",
                    "empty": "0", "invalid": "!"}.get(r.status.value, "?")
            logger.info(f"  {icon} {fname}")
            if r.status in (FileStatus.INVALID, FileStatus.EMPTY):
                logger.info(f"      Status: {r.status.value}  |  Reason: {r.error_message}")
                continue
            md_str = f"Metadata: {r.metadata_extracted}/{r.total_records}"
            asm_str = f"Assembly: {r.assembly_extracted}"
            if r.assembly_no_data > 0:
                asm_str += f" (no_data:{r.assembly_no_data})"
            if r.assembly_failed > 0:
                asm_str += f" (failed:{r.assembly_failed})"
            logger.info(f"      Records: {r.total_records}  |  {md_str}  |  {asm_str}  |  {r.elapsed_time:.1f}s")

            if r.metadata_failed > 0 or r.assembly_failed > 0:
                failed_accs = [t.accession for t in r.record_tracks if not t.metadata_ok]
                if failed_accs:
                    shown = failed_accs[:5]
                    more = f" ... +{len(failed_accs)-5}" if len(failed_accs) > 5 else ""
                    logger.info(f"      Failed ACCESSIONs: {', '.join(shown)}{more}")

    def save_report(self, output_dir: str):
        report_path = os.path.join(output_dir, "extraction_report.json")
        report = {
            "summary": {
                "total_files_found": self.total_files_found,
                "gb_files_found": self.gb_files_found,
                "non_gb_files_skipped": self.non_gb_files_skipped,
                "files_succeeded": self.files_succeeded,
                "files_partial": self.files_partial,
                "files_failed": self.files_failed,
                "files_empty": self.files_empty,
                "files_invalid": self.files_invalid,
                "total_records_in_files": self.total_records_in_files,
                "total_metadata_extracted": self.total_metadata_extracted,
                "total_metadata_failed": self.total_metadata_failed,
                "total_assembly_extracted": self.total_assembly_extracted,
                "total_assembly_no_data": self.total_assembly_no_data,
                "total_assembly_failed": self.total_assembly_failed,
                "total_elapsed_time": round(self.total_elapsed_time, 2),
            },
            "files": []
        }
        for r in self.file_results:
            fr = {
                "file_path": r.file_path,
                "file_name": os.path.basename(r.file_path),
                "status": r.status.value,
                "total_records": r.total_records,
                "metadata_extracted": r.metadata_extracted,
                "metadata_failed": r.metadata_failed,
                "assembly_extracted": r.assembly_extracted,
                "assembly_no_data": r.assembly_no_data,
                "assembly_failed": r.assembly_failed,
                "error_message": r.error_message,
                "elapsed_time": round(r.elapsed_time, 2),
                "record_errors": [
                    {"record_index": e.record_index, "accession": e.accession,
                     "error_type": e.error_type, "error_message": e.error_message,
                     "phase": e.phase}
                    for e in r.record_errors
                ],
            }
            report["files"].append(fr)
        try:
            with open(report_path, 'w', encoding='utf-8') as f:
                json.dump(report, f, ensure_ascii=False, indent=2)
            logger.info(f"Report: {report_path}")
        except Exception as e:
            logger.error(f"Failed to save report: {e}")

    def save_error_csv(self, output_dir: str):
        all_errors = self.get_all_record_errors()
        if not all_errors:
            return
        path = os.path.join(output_dir, "record_errors.csv")
        try:
            with open(path, 'w', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                w.writerow(['file_name', 'record_index', 'accession',
                            'error_type', 'error_message', 'phase'])
                for e in all_errors:
                    w.writerow([os.path.basename(e.file_path), e.record_index,
                                e.accession, e.error_type, e.error_message, e.phase])
            logger.info(f"Error CSV: {path} ({len(all_errors)} records)")
        except Exception as e:
            logger.error(f"Failed to save error CSV: {e}")

    def save_tracking_csv(self, output_dir: str):
        all_tracks = []
        for r in self.file_results:
            all_tracks.extend(r.record_tracks)
        if not all_tracks:
            return
        path = os.path.join(output_dir, "record_tracking.csv")
        try:
            with open(path, 'w', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                w.writerow(['file_name', 'record_index', 'accession',
                            'metadata_ok', 'assembly_ok', 'assembly_source',
                            'has_assembly_data', 'error_message'])
                for r in self.file_results:
                    fname = os.path.basename(r.file_path)
                    for t in r.record_tracks:
                        w.writerow([fname, t.record_index, t.accession,
                                    t.metadata_ok, t.assembly_ok,
                                    t.assembly_source, t.has_assembly_data,
                                    t.error_message])
            logger.info(f"Tracking CSV: {path} ({len(all_tracks)} records)")
        except Exception as e:
            logger.error(f"Failed to save tracking CSV: {e}")


# ===============================
# File discovery and validation
# ===============================
def discover_gb_files(input_paths: List[str], recursive: bool = True) -> Tuple[List[str], int, int]:
    gb_files = []
    total_files = 0
    skipped = 0
    for path in input_paths:
        if os.path.isdir(path):
            walker = os.walk(path) if recursive else [(path, [], os.listdir(path))]
            for root, dirs, files in walker:
                for fn in sorted(files):
                    total_files += 1
                    if fn.lower().endswith(GB_EXTENSIONS):
                        gb_files.append(os.path.join(root, fn))
                    else:
                        skipped += 1
        elif os.path.isfile(path):
            total_files += 1
            if path.lower().endswith(GB_EXTENSIONS):
                gb_files.append(path)
            else:
                skipped += 1
        else:
            logger.warning(f"Path not found: {path}")
    logger.info(f"File discovery: {total_files} files, {len(gb_files)} GenBank, {skipped} skipped")
    return gb_files, total_files, skipped


# P1-11: Combined validation + record counting in single file read
def validate_and_estimate_records(file_path: str) -> Tuple[bool, str, int]:
    """Single-pass: validate file header AND count '//' delimiters."""
    try:
        if not os.path.exists(file_path):
            return False, "File does not exist", 0
        sz = os.path.getsize(file_path)
        if sz == 0:
            return False, "Empty file (0 bytes)", 0

        count = 0
        with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
            header = f.read(min(2048, sz))
            if not header.strip():
                return False, "File content empty", 0
            if not header.lstrip().startswith("LOCUS"):
                return False, f"Does not start with LOCUS", 0
            # Continue reading for record count
            f.seek(0)
            for line in f:
                if line.strip() == "//":
                    count += 1
        return True, "", count
    except PermissionError:
        return False, "No read permission", 0
    except Exception as e:
        return False, f"Validation error: {e}", 0


def format_size(size_bytes) -> str:
    s = float(size_bytes) if size_bytes else 0
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if s < 1024:
            return f"{s:.1f} {unit}"
        s /= 1024
    return f"{s:.1f} PB"


def copy_failed_files(failed_files: List[str], output_dir: str):
    if not failed_files:
        return
    failed_dir = os.path.join(output_dir, "failed_extract")
    os.makedirs(failed_dir, exist_ok=True)
    copied = 0
    for fp in failed_files:
        try:
            dest = os.path.join(failed_dir, os.path.basename(fp))
            if os.path.exists(dest):
                base, ext = os.path.splitext(os.path.basename(fp))
                c = 1
                while os.path.exists(dest):
                    dest = os.path.join(failed_dir, f"{base}_{c}{ext}")
                    c += 1
            shutil.copy2(fp, dest)
            copied += 1
        except Exception as e:
            logger.error(f"Failed to copy {fp}: {e}")
    logger.info(f"Archived {copied}/{len(failed_files)} failed files -> {failed_dir}")


def get_memory_mb() -> Optional[float]:
    if psutil:
        try:
            return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
        except Exception:
            pass
    return None


# ===============================
# NCBI taxonomy
# ===============================
def update_ncbi_taxonomy_database():
    global ete3_available, NCBITaxa
    if not ete3_available:
        logger.warning("ete3 not installed")
        return False
    try:
        if NCBITaxa is None:
            from ete3 import NCBITaxa as C
            NCBITaxa = C
        ncbi = NCBITaxa()
        logger.info("Updating NCBI taxonomy database...")
        ncbi.update_taxonomy_database()
        logger.info("Update complete")
        return True
    except Exception as e:
        logger.error(f"Update error: {e}")
        return False


def initialize_ncbi_taxa(interactive: bool = True) -> bool:
    global ete3_available, NCBITaxa
    ete3_available = importlib.util.find_spec("ete3") is not None
    if ete3_available:
        try:
            from ete3 import NCBITaxa
            logger.info("ete3 module loaded")
            return True
        except Exception as e:
            logger.warning(f"ete3 init error: {e}")
            ete3_available = False
    else:
        logger.warning("ete3 not detected, taxonomy unavailable")
        logger.info("Install: pip install ete3")
        if interactive and not sys.flags.interactive:
            try:
                if input("Continue without taxonomy? (y/n, default y): ").strip().lower()[:1] == 'n':
                    sys.exit(0)
            except (EOFError, IndexError):
                pass
    return False


def get_taxonomy_lineage(tax_id_str: str) -> Dict[str, str]:
    global TAXONOMY_CACHE, ete3_available, NCBITaxa, TAXONOMY_ERROR_REPORTED
    if tax_id_str in TAXONOMY_CACHE:
        return TAXONOMY_CACHE[tax_id_str]
    if not ete3_available or not NCBITaxa:
        return {}
    try:
        tax_id = int(tax_id_str)
    except ValueError:
        TAXONOMY_CACHE[tax_id_str] = {}
        return {}
    try:
        ncbi = NCBITaxa()
        lineage = ncbi.get_lineage(tax_id)
        if not lineage:
            TAXONOMY_CACHE[tax_id_str] = {}
            return {}
        rank_dict = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        result = {
            rank_dict[tid].capitalize(): names[tid]
            for tid in lineage
            if rank_dict.get(tid) in ("class", "order", "family", "genus")
        }
    except Exception as e:
        if tax_id_str not in TAXONOMY_ERROR_REPORTED:
            logger.warning(f"TaxonID '{tax_id_str}' query error: {e}")
            TAXONOMY_ERROR_REPORTED.add(tax_id_str)
        result = {}
    TAXONOMY_CACHE[tax_id_str] = result
    return result


# ===============================
# Metadata extraction (Biopython)
# ===============================
def parse_locus(record: SeqRecord) -> Dict[str, str]:
    return {
        "LocusID": record.id,
        "Length": f"{len(record.seq)} bp",
        "MoleculeType": record.annotations.get("molecule_type", "").title(),
        "Topology": record.annotations.get("topology", "").capitalize(),
        "Division": record.annotations.get("data_file_division", ""),
        "Date": record.annotations.get("date", "")
    }


def parse_source_feature(feature: SeqFeature) -> Dict[str, str]:
    data = {}
    for q, vals in feature.qualifiers.items():
        cleaned = [' '.join(v.replace('\n', ' ').split()) for v in vals]
        if q == "db_xref":
            tids = [x.split(":")[1] for x in cleaned if x.startswith("taxon:")]
            if tids:
                data["TaxonID"] = tids[0]
            continue
        data[q] = "; ".join(cleaned)
    return data


def parse_references(refs: list) -> Dict[str, str]:
    d = {}
    for i, ref in enumerate(refs, 1):
        p = f"Ref{i}"
        d[f"{p}Authors"] = ref.authors or ""
        d[f"{p}Title"] = ref.title or ""
        d[f"{p}Journal"] = ref.journal or ""
        d[f"{p}PubMed"] = ref.pubmed_id or ""
        d[f"{p}Remark"] = ref.comment or ""
    return d


def extract_metadata(record: SeqRecord, include_taxonomy: bool = True) -> Dict[str, str]:
    metadata = {}
    metadata.update(parse_locus(record))
    metadata["ACCESSION"] = record.annotations.get("accessions", [""])[0]
    metadata["Version"] = record.annotations.get("sequence_version", "")
    metadata["Definition"] = record.description
    metadata["Keywords"] = "; ".join(record.annotations.get("keywords", []))
    metadata["Organism"] = record.annotations.get("organism", "")
    tax = record.annotations.get("taxonomy", [])
    if tax:
        metadata["Taxonomy"] = "; ".join(tax)

    for feat in record.features:
        if feat.type == "source":
            metadata.update(parse_source_feature(feat))
            break

    if include_taxonomy:
        tid = metadata.get("TaxonID")
        info = get_taxonomy_lineage(tid) if tid else {}
        for rank in ("Class", "Order", "Family", "Genus"):
            metadata[rank] = info.get(rank, "")

    metadata.update(parse_references(record.annotations.get("references", [])))
    return metadata


# ===============================
# Assembly extraction: 3 strategies
# ===============================
ASSEMBLY_FIELDS = ("Assembly Method", "Assembly Name", "Sequencing Technology",
                   "Genome Coverage", "Expected Final Version")


def _clean_value(v) -> str:
    if isinstance(v, list):
        v = "; ".join(str(x) for x in v)
    return ' '.join(str(v).replace('\n', ' ').split()).strip()


def extract_assembly_from_structured_comment(record: SeqRecord) -> Tuple[Optional[Dict[str, str]], str]:
    sc = record.annotations.get("structured_comment", {})
    asm_data = None
    for key in ("Assembly-Data", "Assembly-data", "assembly-data",
                "Genome-Assembly-Data", "Genome-assembly-data"):
        if key in sc:
            asm_data = sc[key]
            break
    if not asm_data or not isinstance(asm_data, dict):
        return None, ""

    result = {}
    key_map = {
        "Assembly Method": ["Assembly Method", "Assembly_Method", "AssemblyMethod"],
        "Assembly Name": ["Assembly Name", "Assembly_Name", "AssemblyName"],
        "Sequencing Technology": ["Sequencing Technology", "Sequencing_Technology",
                                  "SequencingTechnology"],
        "Genome Coverage": ["Genome Coverage", "Genome_Coverage", "GenomeCoverage"],
        "Expected Final Version": ["Expected Final Version", "Expected_Final_Version"],
    }
    for std_name, variants in key_map.items():
        for vk in variants:
            if vk in asm_data:
                val = _clean_value(asm_data[vk])
                if val:
                    result[std_name] = val
                break
    mapped_keys = set()
    for variants in key_map.values():
        mapped_keys.update(variants)
    for k, v in asm_data.items():
        if k not in mapped_keys:
            val = _clean_value(v)
            if val:
                result[k] = val

    if result:
        return result, "structured_comment"
    return None, ""


def extract_assembly_from_comment_regex(record: SeqRecord) -> Tuple[Optional[Dict[str, str]], str]:
    comment = record.annotations.get("comment", "")
    if not comment:
        return None, ""
    result = {}
    patterns = {
        "Assembly Method": _FLEX_AM,
        "Assembly Name": _FLEX_AN,
        "Sequencing Technology": _FLEX_ST,
        "Genome Coverage": _FLEX_GC,
        "Expected Final Version": _FLEX_EL,
    }
    for field_name, pat in patterns.items():
        m = pat.search(comment)
        if m:
            val = m.group(1).strip().rstrip(';').strip()
            if val:
                result[field_name] = val
    if result:
        return result, "comment_text"
    return None, ""


def extract_assembly_from_raw_block(record_text: str) -> Tuple[Optional[Dict[str, str]], str]:
    block = _BLOCK_RE.search(record_text)
    if not block:
        return None, ""
    content = block.group(1)
    result = {}
    patterns = {
        "Assembly Method": _FLEX_AM,
        "Assembly Name": _FLEX_AN,
        "Sequencing Technology": _FLEX_ST,
        "Genome Coverage": _FLEX_GC,
        "Expected Final Version": _FLEX_EL,
    }
    for field_name, pat in patterns.items():
        m = pat.search(content)
        if m:
            val = m.group(1).strip().rstrip(';').strip()
            if val:
                result[field_name] = val
    if result:
        return result, "regex_block"
    return None, ""


def extract_assembly_from_record(record: SeqRecord) -> Tuple[Optional[Dict[str, str]], str]:
    asm, src = extract_assembly_from_structured_comment(record)
    if asm:
        return asm, src
    asm, src = extract_assembly_from_comment_regex(record)
    if asm:
        return asm, src
    comment = record.annotations.get("comment", "")
    if comment:
        asm, src = extract_assembly_from_raw_block(comment)
        if asm:
            return asm, src
    return None, ""


# ===============================
# Column ordering
# ===============================
BASE_METADATA_COLUMNS = [
    "LocusID", "ACCESSION", "Version", "Length", "MoleculeType", "Topology",
    "Division", "Date", "Definition", "Keywords", "Organism", "Taxonomy",
    "organism", "mol_type", "isolate", "host", "db_xref", "country",
    "collection_date", "collected_by", "identified_by", "note",
    "TaxonID", "Class", "Order", "Family", "Genus",
]

ASSEMBLY_CSV_COLUMNS = [
    'ACCESSION', 'Assembly Method', 'Assembly Name', 'Sequencing Technology',
    'Genome Coverage', 'Expected Final Version',
]


def get_ordered_fieldnames(row: Dict, base: List[str]) -> List[str]:
    ordered = [c for c in base if c in row]
    extra = [c for c in row if c not in ordered]
    return ordered + extra


# ===============================
# Unified record processing
# ===============================
def process_record_unified(
    record: SeqRecord,
    include_taxonomy: bool = True,
    do_metadata: bool = True,
    do_assembly: bool = True,
) -> Tuple[Optional[Dict[str, str]], Optional[Dict[str, str]], str, bool]:
    metadata = None
    assembly_data = None
    assembly_source = ""
    has_asm = False

    if do_metadata:
        metadata = extract_metadata(record, include_taxonomy)

    if do_assembly:
        assembly_data, assembly_source = extract_assembly_from_record(record)
        has_asm = assembly_data is not None
        if assembly_data is not None:
            if "ACCESSION" not in assembly_data:
                acc = record.annotations.get("accessions", [""])[0] or record.id
                assembly_data["ACCESSION"] = acc

    return metadata, assembly_data, assembly_source, has_asm


# ===============================
# P1-7: Shared record processing loop
# ===============================
RecordCallback = Callable[[Optional[Dict[str, str]], Optional[Dict[str, str]],
                           str, bool, RecordTrack, 'FileResult'], None]


def _process_record_loop(
    file_path: str,
    callback: RecordCallback,
    include_taxonomy: bool = True,
    do_metadata: bool = True,
    do_assembly: bool = True,
    log_interval: int = 5000,
) -> FileResult:
    """Shared record processing loop for a single GB file.

    callback receives (metadata, assembly_data, assembly_source, has_asm, track, result) for each record.

    Optimized: Single-pass file reading (no pre-validation).
    """
    result = FileResult(file_path=file_path)
    start = time.time()

    # Quick file existence check
    if not os.path.exists(file_path):
        result.status = FileStatus.INVALID
        result.error_message = "File does not exist"
        result.elapsed_time = time.time() - start
        return result

    sz = os.path.getsize(file_path)
    if sz == 0:
        result.status = FileStatus.EMPTY
        result.error_message = "Empty file (0 bytes)"
        result.elapsed_time = time.time() - start
        return result

    logger.info(f"  Size: {format_size(sz)}")

    record_idx = 0
    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
            for record in SeqIO.parse(handle, "genbank"):
                record_idx += 1
                accession = ""
                track = RecordTrack(accession="", record_index=record_idx)

                try:
                    accession = record.annotations.get("accessions", [""])[0] or record.id
                    track.accession = accession

                    metadata, assembly_data, assembly_source, has_asm = process_record_unified(
                        record, include_taxonomy, do_metadata, do_assembly
                    )

                    callback(metadata, assembly_data, assembly_source, has_asm, track, result)

                    if do_metadata and metadata is not None:
                        result.metadata_extracted += 1
                        track.metadata_ok = True
                    if do_assembly:
                        track.has_assembly_data = has_asm
                        if assembly_data is not None:
                            result.assembly_extracted += 1
                            track.assembly_ok = True
                            track.assembly_source = assembly_source
                        elif not has_asm:
                            result.assembly_no_data += 1

                except Exception as e:
                    track.error_message = str(e)
                    if do_metadata and not track.metadata_ok:
                        result.metadata_failed += 1
                        result.add_record_error(RecordError(
                            file_path=file_path, record_index=record_idx,
                            accession=accession, error_type=type(e).__name__,
                            error_message=str(e), phase="metadata"
                        ))
                    if do_assembly and track.has_assembly_data and not track.assembly_ok:
                        result.assembly_failed += 1

                result.record_tracks.append(track)

                # Progress logging
                done = result.metadata_extracted + result.metadata_failed
                if done > 0 and done % log_interval == 0:
                    elapsed = time.time() - start
                    speed = done / elapsed if elapsed > 0 else 0
                    mem = get_memory_mb()
                    mem_s = f", mem:{mem:.0f}MB" if mem else ""
                    logger.info(
                        f"  Progress: {done} ({speed:.0f}/s{mem_s})")

    except Exception as e:
        result.error_message = f"Biopython parse error: {e}"
        logger.error(f"  File-level error: {e}")

    # Set total records after processing
    result.total_records = record_idx

    # Status determination
    if do_metadata:
        if result.metadata_extracted == 0:
            result.status = FileStatus.FAILED
            if not result.error_message:
                result.error_message = "No metadata extracted"
        elif result.metadata_failed > 0:
            result.status = FileStatus.PARTIAL
        else:
            result.status = FileStatus.SUCCESS
    else:
        result.status = FileStatus.SUCCESS

    result.elapsed_time = time.time() - start
    return result


# ===============================
# Stream mode wrapper
# ===============================
def stream_process_single_file(
    file_path: str,
    metadata_csv_path: Optional[str],
    assembly_csv_path: Optional[str],
    include_taxonomy: bool = True,
    do_extract_metadata: bool = True,  # P1-20: renamed from extract_metadata
    do_extract_assembly: bool = True,
    metadata_delimiter: str = ',',
    assembly_delimiter: str = '\t',
    log_interval: int = 5000,
) -> FileResult:
    metadata_fh = None
    metadata_writer = None
    metadata_header_written = False
    assembly_fh = None
    assembly_writer = None

    try:
        if do_extract_metadata and metadata_csv_path:
            metadata_fh = open(metadata_csv_path, 'w', newline='', encoding='utf-8')
        if do_extract_assembly and assembly_csv_path:
            assembly_fh = open(assembly_csv_path, 'w', newline='', encoding='utf-8')
            assembly_writer = csv.DictWriter(
                assembly_fh, fieldnames=ASSEMBLY_CSV_COLUMNS,
                delimiter=assembly_delimiter, extrasaction='ignore'
            )
            assembly_writer.writeheader()

        def stream_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
            nonlocal metadata_header_written, metadata_writer
            # Write metadata
            if do_extract_metadata and metadata is not None and metadata_fh:
                if not metadata_header_written:
                    fns = get_ordered_fieldnames(metadata, BASE_METADATA_COLUMNS)
                    metadata_writer = csv.DictWriter(
                        metadata_fh, fieldnames=fns,
                        delimiter=metadata_delimiter, extrasaction='ignore'
                    )
                    metadata_writer.writeheader()
                    metadata_header_written = True
                for fn in metadata_writer.fieldnames:
                    if fn not in metadata:
                        metadata[fn] = ""
                metadata_writer.writerow(metadata)

            # Write assembly
            if do_extract_assembly and assembly_data is not None and assembly_fh:
                for col in ASSEMBLY_CSV_COLUMNS:
                    if col not in assembly_data:
                        assembly_data[col] = ""
                assembly_writer.writerow(assembly_data)

        result = _process_record_loop(
            file_path, stream_callback, include_taxonomy,
            do_extract_metadata, do_extract_assembly, log_interval
        )
        return result

    finally:
        for fh in (metadata_fh, assembly_fh):
            if fh:
                try:
                    fh.flush()
                    fh.close()
                except Exception:
                    pass


# ===============================
# Memory mode wrapper
# ===============================
def memory_process_single_file(
    file_path: str,
    include_taxonomy: bool = True,
    do_extract_assembly: bool = True,
    do_extract_metadata: bool = True,  # P1-20: renamed
) -> Tuple[FileResult, List[Dict], List[Dict]]:
    md_rows = []
    asm_rows = []

    def memory_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
        if do_extract_metadata and metadata:
            md_rows.append(metadata)
        if do_extract_assembly and assembly_data:
            asm_rows.append(assembly_data)

    result = _process_record_loop(
        file_path, memory_callback, include_taxonomy,
        do_extract_metadata, do_extract_assembly
    )
    return result, md_rows, asm_rows


# ===============================
# Merge metadata + assembly
# ===============================
def combine_and_save_final_csv(md_file: str, asm_file: str, out_file: str, delimiter: str = ','):
    start = time.time()
    try:
        md_df = pd.read_csv(md_file, low_memory=False)
    except Exception as e:
        logger.error(f"Failed to read metadata: {e}")
        return

    asm_df = None
    try:
        asm_df = pd.read_csv(asm_file, sep='\t')
    except (FileNotFoundError, pd.errors.EmptyDataError):
        pass
    except Exception as e:
        logger.warning(f"Failed to read assembly: {e}")

    if asm_df is not None and "ACCESSION" in md_df.columns and "ACCESSION" in asm_df.columns:
        md_df["ACCESSION"] = md_df["ACCESSION"].astype(str).str.strip()
        asm_df["ACCESSION"] = asm_df["ACCESSION"].astype(str).str.strip()
        final = pd.merge(md_df, asm_df, how="left", on="ACCESSION")
    else:
        final = md_df

    final.to_csv(out_file, index=False, sep=delimiter, encoding='utf-8')
    logger.info(f"Merged: {out_file} ({len(final)} records, {time.time()-start:.1f}s)")


# ===============================
# Main flow: non-batch
# ===============================
def process_all_files(
    gb_files: List[str], output_dir: str,
    metadata_output: str = "metadata.csv",
    assembly_output: str = "assembly.csv",
    final_output: str = "final.csv",
    include_taxonomy: bool = True,
    only_assembly: bool = False,
    only_except_assembly: bool = False,
    metadata_delimiter: str = ',',
    assembly_delimiter: str = '\t',
    final_delimiter: str = ',',
    stream_mode: bool = False,
    log_interval: int = 5000,
) -> ProcessingSummary:
    summary = ProcessingSummary()
    summary.gb_files_found = len(gb_files)
    t0 = time.time()

    do_metadata = not only_assembly
    do_assembly = not only_except_assembly
    tax = include_taxonomy and not only_except_assembly

    metadata_path = os.path.join(output_dir, metadata_output)
    assembly_path = os.path.join(output_dir, assembly_output)
    final_path = os.path.join(output_dir, final_output)

    if stream_mode:
        metadata_fh = None
        metadata_writer = None
        metadata_header_written = False
        assembly_fh = None
        assembly_writer = None
        try:
            if do_metadata:
                metadata_fh = open(metadata_path, 'w', newline='', encoding='utf-8')
            if do_assembly:
                assembly_fh = open(assembly_path, 'w', newline='', encoding='utf-8')
                assembly_writer = csv.DictWriter(
                    assembly_fh, fieldnames=ASSEMBLY_CSV_COLUMNS,
                    delimiter=assembly_delimiter, extrasaction='ignore'
                )
                assembly_writer.writeheader()

            for i, fp in enumerate(gb_files, 1):
                logger.info(f"[{i}/{len(gb_files)}] {os.path.basename(fp)} ({format_size(os.path.getsize(fp))})")

                def stream_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
                    nonlocal metadata_header_written, metadata_writer
                    if do_metadata and metadata and metadata_fh:
                        if not metadata_header_written:
                            fns = get_ordered_fieldnames(metadata, BASE_METADATA_COLUMNS)
                            metadata_writer = csv.DictWriter(
                                metadata_fh, fieldnames=fns,
                                delimiter=metadata_delimiter, extrasaction='ignore')
                            metadata_writer.writeheader()
                            metadata_header_written = True
                        for fn in metadata_writer.fieldnames:
                            if fn not in metadata:
                                metadata[fn] = ""
                        metadata_writer.writerow(metadata)
                    if do_assembly and assembly_data and assembly_fh:
                        for c in ASSEMBLY_CSV_COLUMNS:
                            if c not in assembly_data:
                                assembly_data[c] = ""
                        assembly_writer.writerow(assembly_data)

                fr = _process_record_loop(fp, stream_callback, tax, do_metadata, do_assembly, log_interval)
                summary.add_result(fr)

                icon = {"success": "OK", "partial": "~", "failed": "X"}.get(fr.status.value, "?")
                logger.info(f"  {icon} MD:{fr.metadata_extracted}/{fr.total_records} "
                            f"ASM:{fr.assembly_extracted} ({fr.elapsed_time:.1f}s)")

        finally:
            for fh in (metadata_fh, assembly_fh):
                if fh:
                    try:
                        fh.flush()
                        fh.close()
                    except Exception:
                        pass

        if do_metadata and do_assembly and summary.total_metadata_extracted > 0:
            try:
                combine_and_save_final_csv(metadata_path, assembly_path, final_path, final_delimiter)
            except Exception as e:
                logger.error(f"Merge error: {e}")

    else:
        all_metadata = []
        all_assembly = []
        for i, fp in enumerate(gb_files, 1):
            logger.info(f"[{i}/{len(gb_files)}] {os.path.basename(fp)}")

            def memory_callback(metadata, assembly_data, assembly_source, has_asm, track, result):
                if do_metadata and metadata:
                    all_metadata.append(metadata)
                if do_assembly and assembly_data:
                    all_assembly.append(assembly_data)

            fr = _process_record_loop(fp, memory_callback, tax, do_metadata, do_assembly, log_interval)
            summary.add_result(fr)

            icon = {"success": "OK", "partial": "~", "failed": "X"}.get(fr.status.value, "?")
            logger.info(f"  {icon} MD:{fr.metadata_extracted}/{fr.total_records} "
                        f"ASM:{fr.assembly_extracted} ({fr.elapsed_time:.1f}s)")

        if do_metadata and all_metadata:
            try:
                pd.DataFrame(all_metadata).to_csv(metadata_path, index=False, sep=metadata_delimiter, encoding='utf-8')
                logger.info(f"Metadata: {metadata_path} ({len(all_metadata)} records)")
            except Exception as e:
                logger.error(f"Failed to save metadata: {e}")

        if do_assembly and all_assembly:
            try:
                df = pd.DataFrame(all_assembly)
                for c in ASSEMBLY_CSV_COLUMNS:
                    if c not in df.columns:
                        df[c] = ""
                df[ASSEMBLY_CSV_COLUMNS].to_csv(assembly_path, index=False, sep=assembly_delimiter, encoding='utf-8')
                logger.info(f"Assembly: {assembly_path} ({len(all_assembly)} records)")
            except Exception as e:
                logger.error(f"Failed to save assembly: {e}")

        if do_metadata and do_assembly and all_metadata:
            try:
                combine_and_save_final_csv(metadata_path, assembly_path, final_path, final_delimiter)
            except Exception as e:
                logger.error(f"Merge error: {e}")

    summary.total_elapsed_time = time.time() - t0
    return summary


# ===============================
# Batch mode
# ===============================
def _init_worker(cache_dict):
    """P1-10: Worker initializer for shared taxonomy cache."""
    global TAXONOMY_CACHE
    TAXONOMY_CACHE = cache_dict


def process_batch_mode(
    gb_files: List[str], output_dir: str,
    max_tasks: int = 4, stream_mode: bool = False,
    only_assembly: bool = False, only_except_assembly: bool = False,
    include_taxonomy: bool = True,
    metadata_delimiter: str = ',', assembly_delimiter: str = '\t',
    final_delimiter: str = ',', log_interval: int = 5000,
) -> ProcessingSummary:
    summary = ProcessingSummary()
    summary.gb_files_found = len(gb_files)
    t0 = time.time()

    do_metadata = not only_assembly
    do_assembly = not only_except_assembly
    tax = include_taxonomy and not only_except_assembly

    def _process_single_file_batch(file_path: str) -> FileResult:
        base = os.path.splitext(os.path.basename(file_path))[0]
        sub = os.path.join(output_dir, base + "_metadata")
        os.makedirs(sub, exist_ok=True)

        metadata_path = os.path.join(sub, base + "_metadata.csv") if do_metadata else None
        assembly_path = os.path.join(sub, base + "_assembly.csv") if do_assembly else None
        final_path = os.path.join(sub, base + "_metadata_final.csv")

        if stream_mode:
            fr = stream_process_single_file(
                file_path, metadata_path, assembly_path, tax, do_metadata, do_assembly,
                metadata_delimiter, assembly_delimiter, log_interval)
        else:
            fr, mdr, asmr = memory_process_single_file(file_path, tax, do_assembly, do_metadata)
            if do_metadata and mdr and metadata_path:
                try:
                    pd.DataFrame(mdr).to_csv(metadata_path, index=False, sep=metadata_delimiter, encoding='utf-8')
                except Exception as e:
                    logger.error(f"Save error: {e}")
            if do_assembly and asmr and assembly_path:
                try:
                    df = pd.DataFrame(asmr)
                    for c in ASSEMBLY_CSV_COLUMNS:
                        if c not in df.columns:
                            df[c] = ""
                    df[ASSEMBLY_CSV_COLUMNS].to_csv(assembly_path, index=False, sep=assembly_delimiter, encoding='utf-8')
                except Exception as e:
                    logger.error(f"Save error: {e}")

        if do_metadata and do_assembly and metadata_path and assembly_path:
            if os.path.exists(metadata_path) and os.path.getsize(metadata_path) > 0:
                try:
                    combine_and_save_final_csv(metadata_path, assembly_path, final_path, final_delimiter)
                except Exception as e:
                    logger.error(f"Merge error: {e}")

        return fr

    # P1-10: Shared taxonomy cache via multiprocessing.Manager
    manager = multiprocessing.Manager()
    shared_cache = manager.dict()
    # Pre-populate with current cache
    shared_cache.update(TAXONOMY_CACHE)

    logger.info(f"Batch mode: {len(gb_files)} files, concurrency {max_tasks}, "
                f"stream={'yes' if stream_mode else 'no'}")

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=max_tasks,
        initializer=_init_worker,
        initargs=(shared_cache,)
    ) as ex:
        future_to_file_map = {ex.submit(_process_single_file_batch, f): f for f in gb_files}
        done = 0
        for fut in concurrent.futures.as_completed(future_to_file_map):
            f = future_to_file_map[fut]
            done += 1
            try:
                fr = fut.result()
                summary.add_result(fr)
                icon = {"success": "OK", "partial": "~", "failed": "X"}.get(fr.status.value, "?")
                logger.info(f"[{done}/{len(gb_files)}] {icon} {os.path.basename(f)} "
                            f"MD:{fr.metadata_extracted}/{fr.total_records} "
                            f"ASM:{fr.assembly_extracted} ({fr.elapsed_time:.1f}s)")
            except Exception as e:
                fr = FileResult(file_path=f, status=FileStatus.FAILED, error_message=str(e))
                summary.add_result(fr)
                logger.error(f"[{done}/{len(gb_files)}] X {os.path.basename(f)}: {e}")

    # Sync cache back
    TAXONOMY_CACHE.update(shared_cache)

    summary.total_elapsed_time = time.time() - t0
    return summary


# ===============================
# Public API
# ===============================
def extract(
    input_files: list,
    output_dir: str,
    stream: bool = True,
    batch: bool = False,
    max_tasks: int = 1,
    final_name: str = "final.csv",
) -> StepResult:
    """Extract metadata from GenBank files. Returns StepResult."""
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    try:
        gb_files, total, skipped = discover_gb_files(input_files, True)
        if not gb_files:
            return StepResult(success=False, output_file="", elapsed=time.time() - start_time)

        valid_files = []
        for f in gb_files:
            ok, err, _ = validate_and_estimate_records(f)
            if ok:
                valid_files.append(f)

        if not valid_files:
            return StepResult(success=False, output_file="", elapsed=time.time() - start_time)

        initialize_ncbi_taxa()

        if batch:
            summary = process_batch_mode(
                valid_files, output_dir, max_tasks, stream,
                False, False, True, ',', '\t', ',', 5000)
        else:
            summary = process_all_files(
                valid_files, output_dir,
                "metadata.csv", "assembly.csv", final_name,
                True, False, False, ',', '\t', ',', stream, 5000)

        final_path = os.path.join(output_dir, final_name)
        rows = 0
        if os.path.exists(final_path):
            try:
                rows = len(pd.read_csv(final_path, dtype=str))
            except Exception:
                pass

        return StepResult(
            success=True,
            output_file=final_path,
            rows=rows,
            elapsed=time.time() - start_time,
        )
    except Exception as e:
        logger.error(f"Extract failed: {e}")
        return StepResult(success=False, output_file="", elapsed=time.time() - start_time)


# ===============================
# Main function
# ===============================
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="GenBank Metadata Extractor (refactored): 3-strategy Assembly + record tracking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, nargs="+",
                        help="Input GenBank files or directories")
    parser.add_argument("-o", "--output", default=".", help="Output directory")
    parser.add_argument("-m", "--metadata-output", default="metadata.csv")
    parser.add_argument("-a", "--assembly-output", default="assembly.csv")
    parser.add_argument("-f", "--final-output", default="final.csv")

    parser.add_argument("--stream", action="store_true",
                        help="Stream processing (recommended for large files)")
    parser.add_argument("--batch", action="store_true",
                        help="Batch mode (parallel)")
    parser.add_argument("--max_tasks", type=int, default=4)
    parser.add_argument("--log_interval", type=int, default=5000)

    parser.add_argument("--metadata_delimiter", default=",")
    parser.add_argument("--assembly_delimiter", default="\t")
    parser.add_argument("--final_delimiter", default=",")

    exclusive_group = parser.add_mutually_exclusive_group()
    exclusive_group.add_argument("--only_assembly_meta", "-oam", action="store_true")
    exclusive_group.add_argument("--only_except_assembly_meta", "-oeam", action="store_true")

    parser.add_argument("--update_taxonomy", action="store_true")
    parser.add_argument("--no_report", action="store_true")
    parser.add_argument("--no_tracking", action="store_true")
    parser.add_argument("--no_recursive", action="store_true")
    parser.add_argument("--if_use_chunksize", action="store_true",
                        help="(Legacy) same as --stream")

    args = parser.parse_args(argv)
    stream = args.stream or args.if_use_chunksize

    if psutil:
        mem = psutil.virtual_memory()
        logger.info(f"Available memory: {mem.available / (1024**3):.2f} GB")

    os.makedirs(args.output, exist_ok=True)

    # Step 1: Discover files
    logger.info("=" * 75)
    logger.info("Step 1/4: Discover GenBank files")
    logger.info("=" * 75)
    gb_files, total, skipped = discover_gb_files(args.input, not args.no_recursive)

    if not gb_files:
        logger.error("No GenBank files found")
        sys.exit(1)

    valid_files = []
    summary = ProcessingSummary()
    summary.total_files_found = total
    summary.non_gb_files_skipped = skipped
    summary.gb_files_found = len(gb_files)

    for f in gb_files:
        ok, err, _ = validate_and_estimate_records(f)
        if ok:
            valid_files.append(f)
        else:
            fr = FileResult(file_path=f, status=FileStatus.INVALID, error_message=err)
            summary.add_result(fr)
            logger.warning(f"  Invalid: {os.path.basename(f)} ({err})")

    logger.info(f"Valid files: {len(valid_files)}/{len(gb_files)}")

    # Pre-scan
    total_bytes = sum(os.path.getsize(f) for f in valid_files if os.path.exists(f))
    if not stream and total_bytes > 500 * 1024 * 1024:
        logger.warning(f"Large dataset ({format_size(total_bytes)}), recommend --stream")

    if not valid_files:
        copy_failed_files(summary.get_failed_files(), args.output)
        if not args.no_report:
            summary.save_report(args.output)
        summary.print_summary()
        sys.exit(1)

    # Step 2: Taxonomy
    if not args.only_assembly_meta and not args.only_except_assembly_meta:
        logger.info("=" * 75)
        logger.info("Step 2/4: Taxonomy module")
        logger.info("=" * 75)
        initialize_ncbi_taxa()
        if args.update_taxonomy:
            update_ncbi_taxonomy_database()

    # Step 3: Extract
    logger.info("=" * 75)
    mode_str = f"{'stream' if stream else 'memory'} / {'batch' if args.batch else 'merged'}"
    logger.info(f"Step 3/4: Extract metadata ({mode_str})")
    logger.info("=" * 75)

    if args.batch:
        processing_result = process_batch_mode(
            valid_files, args.output, args.max_tasks, stream,
            args.only_assembly_meta, args.only_except_assembly_meta,
            not args.only_except_assembly_meta,
            args.metadata_delimiter, args.assembly_delimiter,
            args.final_delimiter, args.log_interval)
        for r in processing_result.file_results:
            summary.add_result(r)
        summary.total_elapsed_time = processing_result.total_elapsed_time
    else:
        processing_result = process_all_files(
            valid_files, args.output,
            args.metadata_output, args.assembly_output, args.final_output,
            not args.only_except_assembly_meta,
            args.only_assembly_meta, args.only_except_assembly_meta,
            args.metadata_delimiter, args.assembly_delimiter,
            args.final_delimiter, stream, args.log_interval)
        for r in processing_result.file_results:
            summary.add_result(r)
        summary.total_elapsed_time = processing_result.total_elapsed_time

    # Step 4: Report
    logger.info("=" * 75)
    logger.info("Step 4/4: Report & archive")
    logger.info("=" * 75)

    failed = summary.get_failed_files()
    if failed:
        copy_failed_files(failed, args.output)

    summary.print_summary()

    if not args.no_report:
        summary.save_report(args.output)
        summary.save_error_csv(args.output)
    if not args.no_tracking:
        summary.save_tracking_csv(args.output)

    if summary.files_failed > 0 or summary.files_invalid > 0:
        logger.info("Some files failed, check failed_extract/ and extraction_report.json")
    else:
        logger.info("All processing complete!")


if __name__ == "__main__":
    main()
