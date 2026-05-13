"""Tests for classify.py matching logic — pure functions only."""

import pytest

from g2t.classify import (
    MatchConfig,
    _BUILTIN_GENE_DICT,
    any_match,
    categorize_unmatched_definition,
    clean_length,
    conflict_check_and_reassign,
    first_round_match,
    init_gene_groups,
    match_mito_genes,
    match_nuclear_genes,
    recheck_match_row,
    smart_match,
)


@pytest.fixture
def config():
    return MatchConfig()


# --- smart_match ---


class TestSmartMatch:
    def test_case_insensitive(self):
        # smart_match looks for the phrase (from GENE_DICT) in text; "coi" is a keyword
        assert smart_match("coi", "COI gene sequence") is True
        assert smart_match("cytochrome oxidase subunit I", "CYTOCHROME OXIDASE SUBUNIT I") is True

    def test_word_boundary(self):
        assert smart_match("coi", "coincidence factor") is False

    def test_exact_phrase(self):
        assert smart_match("cytochrome b", "cytochrome b gene") is True

    def test_no_match(self):
        assert smart_match("coi", "18S ribosomal RNA") is False

    @pytest.mark.parametrize("gene", ["coi", "16s", "12s", "cob", "28s", "18s", "h3"])
    def test_each_gene_type_has_matching_synonyms(self, gene):
        """Verify each gene type has at least 3 synonyms that match."""
        keywords = _BUILTIN_GENE_DICT.get(gene, [])
        assert len(keywords) >= 3, f"{gene} has fewer than 3 synonyms"
        matched_count = sum(1 for kw in keywords if smart_match(kw, kw))
        assert matched_count >= 3


# --- any_match ---


class TestAnyMatch:
    def test_match_found(self):
        assert any_match(["coi", "cox1"], "COI gene sequence") is True

    def test_no_match(self):
        assert any_match(["coi", "cox1"], "18S ribosomal RNA") is False

    def test_empty_list(self):
        assert any_match([], "anything") is False


# --- clean_length ---


class TestCleanLength:
    def test_with_bp(self):
        assert clean_length("650 bp") == 650

    def test_without_bp(self):
        assert clean_length("650") == 650

    def test_none(self):
        assert clean_length(None) == 0

    def test_empty(self):
        assert clean_length("") == 0

    def test_string_with_spaces(self):
        assert clean_length(" 1200 bp ") == 1200


# --- match_mito_genes ---


class TestMatchMitoGenes:
    def test_coi_match(self, config):
        matched, notes = match_mito_genes(
            "cytochrome oxidase subunit I partial sequence", 650, "dna", config
        )
        assert "coi" in matched

    def test_mtgenome_priority(self, config):
        matched, notes = match_mito_genes(
            "mitochondrion, complete genome; cytochrome oxidase subunit I", 16000, "dna", config
        )
        assert "mtgenome" in matched
        assert "coi" not in matched  # mtgenome takes priority

    def test_length_out_of_range(self, config):
        matched, notes = match_mito_genes("cytochrome oxidase subunit I", 10, "dna", config)
        assert len(matched) == 0

    def test_16s_match(self, config):
        matched, notes = match_mito_genes(
            "16S ribosomal RNA gene, partial sequence; mitochondrial", 550, "dna", config
        )
        assert "16s" in matched


# --- match_nuclear_genes ---


class TestMatchNuclearGenes:
    def test_18s_match(self, config):
        matched, notes = match_nuclear_genes(
            "18S ribosomal RNA, partial sequence", 1800, config
        )
        assert "18s" in matched

    def test_28s_match(self, config):
        matched, notes = match_nuclear_genes(
            "28S ribosomal RNA gene, partial sequence", 800, config
        )
        assert "28s" in matched

    def test_h3_match(self, config):
        matched, notes = match_nuclear_genes(
            "histone H3 gene, partial cds", 350, config
        )
        assert "h3" in matched

    def test_18_28s_priority(self, config):
        # Use the exact synonym from _BUILTIN_GENE_DICT for 18-28s
        matched, notes = match_nuclear_genes(
            "18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, "
            "5.8S rRNA gene, and internal transcribed spacer 2, complete sequence; "
            "and 28S ribosomal RNA gene, partial sequence",
            8000,  # within 18-28s length range (5000-15000)
            config,
        )
        assert "18-28s" in matched


# --- first_round_match ---


class TestFirstRoundMatch:
    def test_mito_organelle_routes_to_mito(self, config):
        row = {
            "Definition": "cytochrome oxidase subunit I",
            "organelle": "mitochondrion",
            "Length": "650",
            "MoleculeType": "genomic DNA",
            "LocusID": "LOC1",
        }
        result = first_round_match(row, config)
        assert "coi" in result["gene_type"]

    def test_nuclear_organelle_routes_to_nuclear(self, config):
        row = {
            "Definition": "18S ribosomal RNA, partial sequence",
            "organelle": "nucleus",
            "Length": "1800",
            "MoleculeType": "genomic DNA",
            "LocusID": "LOC2",
        }
        result = first_round_match(row, config)
        assert "18s" in result["gene_type"]

    def test_no_match_returns_empty(self, config):
        row = {
            "Definition": "unknown sequence fragment",
            "organelle": "",
            "Length": "100",
            "MoleculeType": "genomic DNA",
            "LocusID": "LOC3",
        }
        result = first_round_match(row, config)
        assert result["gene_type"] == ""


# --- conflict_check_and_reassign ---


class TestConflictCheck:
    def test_mtgenome_skips_conflict(self):
        result = conflict_check_and_reassign(
            "mtgenome", "mitochondrion complete genome", "circular", "mito", 16000
        )
        assert "mtgenome" in result["gene_type"]
        assert not result["conflict"]

    def test_12s_with_large_subunit_conflict(self):
        result = conflict_check_and_reassign(
            "12s", "large subunit ribosomal RNA gene", "linear", "mito", 550
        )
        assert result["conflict"] is True

    def test_16s_with_small_subunit_conflict(self):
        result = conflict_check_and_reassign(
            "16s", "small subunit ribosomal RNA gene", "linear", "mito", 550
        )
        assert result["conflict"] is True

    def test_multi_mito_three_genes_long_reassigns_mtgenome(self):
        # Conflict check: >2 mito genes + Length>3000 -> mtgenome
        result = conflict_check_and_reassign(
            "coi,16s,12s", "some definition", "linear", "mito", 5000
        )
        assert result["gene_type"] == "mtgenome"

    def test_18_28s_skips_conflict(self):
        result = conflict_check_and_reassign(
            "18-28s", "some definition", "linear", "nuclear", 8000
        )
        assert "18-28s" in result["gene_type"]
        assert not result["conflict"]

    def test_single_gene_no_conflict(self):
        result = conflict_check_and_reassign(
            "coi", "cytochrome oxidase subunit I", "linear", "mito", 650
        )
        assert "coi" in result["gene_type"]
        assert not result["conflict"]


# --- recheck_match_row ---


class TestRecheckMatchRow:
    def test_circular_topology_gives_mtgenome(self, config):
        row = {
            "Definition": "mitochondrion complete genome",
            "organelle": "mitochondrion",
            "MoleculeType": "dna",
            "Topology": "circular",
            "Length": "16000",
        }
        result = recheck_match_row(row, config)
        assert "mtgenome" in result["gene_type"]

    def test_pcr_primer_coi(self, config):
        row = {
            "Definition": "unknown gene",
            "organelle": "",
            "MoleculeType": "dna",
            "Topology": "linear",
            "Length": "650",
            "PCR_primers": "fwd_name: LCO1490, rev_name: HCO2198",
        }
        result = recheck_match_row(row, config)
        assert "coi" in result["gene_type"]

    def test_no_match_returns_empty(self, config):
        row = {
            "Definition": "unknown hypothetical protein",
            "organelle": "",
            "MoleculeType": "dna",
            "Topology": "linear",
            "Length": "300",
        }
        result = recheck_match_row(row, config)
        assert result["gene_type"] == ""


# --- categorize_unmatched_definition ---


class TestCategorizeUnmatchedDefinition:
    def test_whole_genome(self):
        assert categorize_unmatched_definition("whole genome shotgun sequence") == "whole_genome_shotgun"

    def test_mrna(self):
        # Checks for 'mrna sequence' or 'cdna clone'
        assert categorize_unmatched_definition("mrna sequence for something") == "mRNA/cDNA"
        assert categorize_unmatched_definition("cDNA clone ABC123") == "mRNA/cDNA"

    def test_nadh(self):
        assert categorize_unmatched_definition("NADH dehydrogenase subunit 5") == "NADH_dehydrogenase"

    def test_atp_synthase(self):
        assert categorize_unmatched_definition("ATP synthase subunit") == "ATP_synthase"

    def test_other(self):
        assert categorize_unmatched_definition("random unknown sequence") == "other"

    def test_trna(self):
        assert categorize_unmatched_definition("tRNA-Leu gene") == "tRNA"

    def test_dead_code_nadh_branch(self):
        """Line 1005 'elif nad' already catches 'nadh', making line 1013 'elif nadh' unreachable."""
        result = categorize_unmatched_definition("nadh dehydrogenase")
        assert result == "NADH_dehydrogenase"
        # Both "nad" and "nadh" in def_lower, but only line 1005 fires
