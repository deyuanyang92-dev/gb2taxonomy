"""Round 1: match a record's DEFINITION to gene types per compartment."""

from __future__ import annotations

import logging
import re
from typing import Any

from g2t._classify.config import MatchConfig
from g2t._classify.text import _negative_hit, any_match, clean_definition, clean_length
from g2t.utils import length_in_range, parse_interval

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


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
