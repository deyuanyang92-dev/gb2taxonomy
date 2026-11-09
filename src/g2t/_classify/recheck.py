"""Round 2: recheck unmatched records against all compartments."""

from __future__ import annotations

import logging
import re
from typing import Any

from g2t._classify.config import MatchConfig
from g2t._classify.match import _build_result, _mtgenome_min, _record_group
from g2t._classify.text import _negative_hit, any_match, clean_definition, clean_length

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


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
