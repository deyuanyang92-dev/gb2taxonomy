"""One GenBank record -> metadata row and assembly row; column order."""

from __future__ import annotations

import logging
import re
from typing import Any

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from g2t._extract.files import open_genbank
from g2t._extract.taxonomy import get_taxonomy_lineage

logger = logging.getLogger("g2t.extract")  # one logger name for the whole step


# ===============================
# Assembly regex (pre-compiled)
# ===============================
# a value ends at the line end or at the next field label on the same line ("...: Illumina; Assembly Method: ...")
_FLEX_END = (
    r"(.*?)(?=;\s*(?:Assembly\s+Method|Assembly\s+Name|Sequencing\s+Technology|Genome\s+Coverage|"
    r"Expected\s+Final\s+Version)\s*::?|\n|$)"
)
_FLEX_AM = re.compile(r"Assembly\s+Method\s*::?\s*" + _FLEX_END, re.IGNORECASE)
_FLEX_AN = re.compile(r"Assembly\s+Name\s*::?\s*" + _FLEX_END, re.IGNORECASE)
_FLEX_ST = re.compile(r"Sequencing\s+Technology\s*::?\s*" + _FLEX_END, re.IGNORECASE)
_FLEX_GC = re.compile(r"Genome\s+Coverage\s*::?\s*" + _FLEX_END, re.IGNORECASE)
_FLEX_EL = re.compile(r"Expected\s+Final\s+Version\s*::?\s*" + _FLEX_END, re.IGNORECASE)


# ===============================
# Metadata extraction (Biopython)
# ===============================
def parse_locus(record: SeqRecord) -> dict[str, str]:
    ann = record.annotations
    return {
        "LocusID": str(record.id),
        "Length": f"{len(record.seq) if record.seq is not None else 0} bp",
        "MoleculeType": str(ann.get("molecule_type", "") or "").title(),
        "Topology": str(ann.get("topology", "") or "").capitalize(),
        "Division": str(ann.get("data_file_division", "") or ""),
        "Date": str(ann.get("date", "") or ""),
    }


def parse_source_feature(feature: SeqFeature) -> dict[str, str]:
    data = {}
    for q, vals in feature.qualifiers.items():
        cleaned = [" ".join(v.replace("\n", " ").split()) for v in vals]
        if q == "db_xref":
            # The first taxon:<id> becomes TaxonID.  Every other cross-reference
            # (BOLD:..., a second taxon:..., etc.) is kept in the db_xref column.
            tids = [x.split(":", 1)[1].strip() for x in cleaned if x.startswith("taxon:")]
            tids = [t for t in tids if t]
            if tids:
                data["TaxonID"] = tids[0]
            others = [x for x in cleaned if not x.startswith("taxon:")]
            others += [f"taxon:{t}" for t in tids[1:]]
            if others:
                data["db_xref"] = "; ".join(others)
            continue
        data[q] = "; ".join(cleaned)
    return data


def parse_references(refs: list) -> dict[str, str]:
    d = {}
    for i, ref in enumerate(refs, 1):
        p = f"Ref{i}"
        d[f"{p}Authors"] = ref.authors or ""
        d[f"{p}Title"] = ref.title or ""
        d[f"{p}Journal"] = ref.journal or ""
        d[f"{p}PubMed"] = ref.pubmed_id or ""
        d[f"{p}Remark"] = ref.comment or ""
    return d


def extract_metadata(record: SeqRecord, include_taxonomy: bool = True) -> dict[str, str]:
    ann: dict[str, Any] = record.annotations
    metadata: dict[str, str] = {}
    metadata.update(parse_locus(record))
    metadata["ACCESSION"] = str((ann.get("accessions") or [""])[0])
    metadata["Version"] = str(ann.get("sequence_version", "") or "")
    metadata["Definition"] = record.description
    metadata["Keywords"] = "; ".join(map(str, ann.get("keywords") or []))
    metadata["Organism"] = str(ann.get("organism", "") or "")
    tax = ann.get("taxonomy") or []
    if tax:
        metadata["Taxonomy"] = "; ".join(map(str, tax))

    sources = [feat for feat in record.features if feat.type == "source"]
    if len(sources) > 1:
        logger.warning(
            f"{metadata['ACCESSION'] or record.id}: {len(sources)} source features found; "
            f"using the first one and ignoring the others"
        )
    if sources:
        metadata.update(parse_source_feature(sources[0]))
    if not metadata.get("organism"):
        # No source feature (or no /organism in it): fall back to the record-level organism
        metadata["organism"] = metadata.get("Organism", "")

    if include_taxonomy:
        tid = metadata.get("TaxonID")
        info = get_taxonomy_lineage(tid) if tid else {}
        for rank in ("Class", "Order", "Family", "Genus"):
            metadata[rank] = info.get(rank, "")

    metadata.update(parse_references(ann.get("references") or []))
    metadata.update(parse_header_extras(record))
    if metadata.get("PCR_primers"):
        metadata.update(parse_primers(metadata["PCR_primers"]))
    return metadata


def parse_primers(raw: str) -> dict[str, str]:
    """Primer names from /PCR_primers ("fwd_name: LCO1490, fwd_seq: ggtc..., rev_name: HCO2198, ...";
    several sets separated by ';'): primer_fwd / primer_rev, names joined with '; '."""
    out = {}
    for side in ("fwd", "rev"):
        names = re.findall(rf"{side}_name:\s*([^,;]+)", raw)
        seen = list(dict.fromkeys(n.strip() for n in names if n.strip()))
        if seen:
            out[f"primer_{side}"] = "; ".join(seen)
    return out


def parse_header_extras(record: SeqRecord) -> dict[str, str]:
    """Header fields beyond the source feature and references: the free-text COMMENT (where submitters
    name the platform, e.g. "based on 91x PacBio data"), DBLINK (BioProject, BioSample, SRA runs) and
    structured comments other than Assembly-Data ("##MIGS-Data-START##" ...)."""
    out: dict[str, str] = {}
    comment = " ".join(str(record.annotations.get("comment", "") or "").split())
    if comment:
        out["Comment"] = comment
    links: dict[str, list[str]] = {}
    for x in record.dbxrefs or []:
        k, _, v = str(x).partition(":")
        key = {"Sequence Read Archive": "SRA"}.get(k.strip(), k.strip())
        if key in ("BioProject", "BioSample", "SRA", "Assembly"):
            links.setdefault(key, []).extend(i.strip() for i in v.split(",") if i.strip())
    for k, ids in links.items():
        out[k] = "; ".join(dict.fromkeys(ids))
    sc: Any = record.annotations.get("structured_comment") or {}
    other = [f"{block}: " + "; ".join(f"{k}={v}" for k, v in fields.items())
             for block, fields in (sc.items() if isinstance(sc, dict) else [])
             if isinstance(fields, dict) and "assembly-data" not in str(block).lower()]
    if other:
        out["StructuredComment"] = " | ".join(other)
    return out


# ===============================
# Assembly extraction: 3 strategies
# ===============================
ASSEMBLY_FIELDS = (
    "Assembly Method",
    "Assembly Name",
    "Sequencing Technology",
    "Genome Coverage",
    "Expected Final Version",
)


def _clean_value(v) -> str:
    if isinstance(v, list):
        v = "; ".join(str(x) for x in v)
    return " ".join(str(v).replace("\n", " ").split()).strip()


def extract_assembly_from_structured_comment(record: SeqRecord) -> tuple[dict[str, str] | None, str]:
    sc: Any = record.annotations.get("structured_comment", {})
    asm_data = None
    for key in ("Assembly-Data", "Assembly-data", "assembly-data", "Genome-Assembly-Data", "Genome-assembly-data"):
        if key in sc:
            asm_data = sc[key]
            break
    if not asm_data or not isinstance(asm_data, dict):
        return None, ""

    result = {}
    key_map = {
        "Assembly Method": ["Assembly Method", "Assembly_Method", "AssemblyMethod"],
        "Assembly Name": ["Assembly Name", "Assembly_Name", "AssemblyName"],
        "Sequencing Technology": ["Sequencing Technology", "Sequencing_Technology", "SequencingTechnology"],
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


def extract_assembly_from_comment_regex(record: SeqRecord) -> tuple[dict[str, str] | None, str]:
    comment = str(record.annotations.get("comment", "") or "")
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
            val = m.group(1).strip().rstrip(";").strip()
            if val:
                result[field_name] = val
    if result:
        return result, "comment_text"
    return None, ""


def extract_assembly_from_record(record: SeqRecord) -> tuple[dict[str, str] | None, str]:
    asm, src = extract_assembly_from_structured_comment(record)
    if asm:
        return asm, src
    asm, src = extract_assembly_from_comment_regex(record)
    if asm:
        return asm, src
    return None, ""


# ===============================
# Column ordering
# ===============================
BASE_METADATA_COLUMNS = [
    "LocusID",
    "ACCESSION",
    "Version",
    "Length",
    "MoleculeType",
    "Topology",
    "Division",
    "Date",
    "Definition",
    "Keywords",
    "Organism",
    "Taxonomy",
    "organism",
    "mol_type",
    "organelle",
    "isolate",
    "host",
    "db_xref",
    "country",
    "specimen_voucher",
    "culture_collection",
    "clone",
    "strain",
    "collection_date",
    "collected_by",
    "identified_by",
    "note",
    "TaxonID",
    "Class",
    "Order",
    "Family",
    "Genus",
]

ASSEMBLY_CSV_COLUMNS = [
    "ACCESSION",
    "Assembly Method",
    "Assembly Name",
    "Sequencing Technology",
    "Genome Coverage",
    "Expected Final Version",
]


def get_ordered_fieldnames(row: dict, base: list[str]) -> list[str]:
    ordered = list(base)
    extra = [c for c in row if c not in ordered]
    return ordered + extra


def canonical_metadata_columns(present) -> list[str]:
    """Column order shared by every output path: BASE columns, then extras sorted.

    Identical to what ``prescan_fieldnames`` returns when ``present`` holds the
    keys that occur in the data, but needs no second parse of the files.
    """
    base = set(BASE_METADATA_COLUMNS)
    return list(BASE_METADATA_COLUMNS) + sorted(c for c in set(present) if c not in base)


def prescan_fieldnames(gb_files: list[str]) -> list[str]:
    """Pre-scan all GB files to collect every unique source qualifier key."""

    all_keys = set()
    for fp in gb_files:
        try:
            with open_genbank(fp) as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    for feat in record.features:
                        if feat.type == "source":
                            for q in feat.qualifiers:
                                if q == "db_xref":
                                    all_keys.add("TaxonID")
                                else:
                                    all_keys.add(q)
                            break
                    all_keys.update(parse_header_extras(record))  # Comment, BioProject, SRA, ...
                    if "PCR_primers" in all_keys:
                        all_keys.update(("primer_fwd", "primer_rev"))
                    # References keys
                    refs = record.annotations.get("references", [])
                    for i in range(1, len(refs) + 1):
                        p = f"Ref{i}"
                        for suffix in ("Authors", "Title", "Journal", "PubMed", "Remark"):
                            all_keys.add(f"{p}{suffix}")
        except Exception as e:
            logger.warning(f"Pre-scan of {fp} stopped early: {e}; some columns may be missing")
    # Build ordered list: BASE first, then discovered extras
    ordered = list(BASE_METADATA_COLUMNS)
    for k in sorted(all_keys):
        if k not in ordered:
            ordered.append(k)
    return ordered
