"""Round 1.5: check the matched gene type against the title and length; reassign."""

from __future__ import annotations

import logging
import re
from typing import Any

from g2t._classify.match import _build_result

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


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
