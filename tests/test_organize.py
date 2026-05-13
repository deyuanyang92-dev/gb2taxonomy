"""Tests for organize.py pure functions and integration."""

import pytest
import pandas as pd
import numpy as np

from g2t.organize import (
    OrganizeConfig,
    build_header,
    check_columns,
    collect_locusids_per_gene,
    get_first_nonempty,
    metadata_advanced,
    metadata_all,
    metadata_first_nonempty,
    metadata_per_gene,
    normalize_gene_name,
    organize,
    process_group,
    _format_gene_values,
    _propagate_locusids,
)


@pytest.fixture
def config():
    return OrganizeConfig()


@pytest.fixture
def sample_group():
    return pd.DataFrame({
        "species_voucher_new": ["sp1", "sp1", "sp1"],
        "organism": ["Homo sapiens", "Homo sapiens", "Homo sapiens"],
        "gene_type": ["coi", "16s", "18s"],
        "LocusID": ["LOC1", "LOC2", "LOC3"],
        "Class": ["Mammalia", "Mammalia", "Mammalia"],
        "Order": ["Primates", "Primates", "Primates"],
    })


# --- OrganizeConfig ---


class TestOrganizeConfig:
    def test_default_gene_order(self):
        cfg = OrganizeConfig()
        assert "coi" in cfg.gene_order
        assert "mtgenome" in cfg.gene_order
        assert len(cfg.gene_order) == 13

    def test_get_mito_genes(self, config):
        mito = config.get_mito_genes()
        assert "coi" in mito
        assert "16s" in mito
        assert "18s" not in mito

    def test_get_mito_genes_excludes_nuclear(self, config):
        mito = config.get_mito_genes()
        for gene in ["18s", "28s", "h3", "ef-1", "its1-its2", "18-28s"]:
            assert gene not in mito

    def test_default_metadata_mode(self):
        cfg = OrganizeConfig()
        assert cfg.metadata_mode == "first_nonempty"


# --- normalize_gene_name ---


class TestNormalizeGeneName:
    def test_plain_name(self):
        assert normalize_gene_name("coi") == "coi"

    def test_trailing_hash(self):
        assert normalize_gene_name("coi#") == "coi"

    def test_trailing_question_mark(self):
        assert normalize_gene_name("16s?") == "16s"

    def test_multiple_trailing(self):
        assert normalize_gene_name("coi##??") == "coi"

    def test_no_change_without_trailing(self):
        assert normalize_gene_name("mtgenome") == "mtgenome"


# --- get_first_nonempty ---


class TestGetFirstNonempty:
    def test_first_value(self):
        s = pd.Series(["a", "b", "c"])
        assert get_first_nonempty(s) == "a"

    def test_skip_nan(self):
        s = pd.Series([np.nan, "b", "c"])
        assert get_first_nonempty(s) == "b"

    def test_skip_empty(self):
        s = pd.Series(["", "  ", "c"])
        assert get_first_nonempty(s) == "c"

    def test_all_empty(self):
        s = pd.Series([np.nan, "", "  "])
        assert get_first_nonempty(s) == ""

    def test_strips_whitespace(self):
        s = pd.Series(["  hello  "])
        assert get_first_nonempty(s) == "hello"


# --- check_columns ---


class TestCheckColumns:
    def test_all_mandatory_present(self):
        df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
        assert check_columns(df, ["a", "b"]) is True

    def test_missing_mandatory(self):
        df = pd.DataFrame({"a": [1]})
        assert check_columns(df, ["a", "b"]) is False

    def test_adds_missing_optional(self):
        df = pd.DataFrame({"a": [1]})
        check_columns(df, ["a"], optional=["b", "c"])
        assert "b" in df.columns
        assert "c" in df.columns
        assert df["b"].iloc[0] == ""

    def test_no_optional_arg(self):
        df = pd.DataFrame({"a": [1]})
        assert check_columns(df, ["a"]) is True


# --- _propagate_locusids ---


class TestPropagateLocusids:
    def test_propagate_mtgenome_to_children(self):
        gene_ids = {
            "mtgenome": ["MG1"],
            "coi": [],
            "16s": [],
            "12s": [],
        }
        _propagate_locusids(gene_ids, "mtgenome", ["coi", "16s", "12s"])
        assert "MG1" in gene_ids["coi"]
        assert "MG1" in gene_ids["16s"]
        assert "MG1" in gene_ids["12s"]

    def test_no_propagation_when_parent_empty(self):
        gene_ids = {"mtgenome": [], "coi": ["C1"]}
        _propagate_locusids(gene_ids, "mtgenome", ["coi"])
        assert gene_ids["coi"] == ["C1"]

    def test_child_preserves_original_order(self):
        gene_ids = {"mtgenome": ["M1", "M2"], "coi": ["C1"]}
        _propagate_locusids(gene_ids, "mtgenome", ["coi"])
        assert gene_ids["coi"] == ["C1", "M1", "M2"]

    def test_no_duplicates(self):
        gene_ids = {"mtgenome": ["X1"], "coi": ["X1"]}
        _propagate_locusids(gene_ids, "mtgenome", ["coi"])
        assert gene_ids["coi"] == ["X1"]


# --- _format_gene_values ---


class TestFormatGeneValues:
    def test_empty(self):
        assert _format_gene_values({}) == ""

    def test_all_empty_sets(self):
        assert _format_gene_values({"coi": set(), "16s": set()}) == ""

    def test_single_value(self):
        assert _format_gene_values({"coi": {"val1"}}) == "val1"

    def test_same_value_multiple_genes(self):
        result = _format_gene_values({"coi": {"val1"}, "16s": {"val1"}})
        assert result == "val1"

    def test_different_values(self):
        result = _format_gene_values({"coi": {"A"}, "16s": {"B"}})
        assert "coi" in result
        assert "16s" in result
        assert "A" in result
        assert "B" in result


# --- collect_locusids_per_gene ---


class TestCollectLocusidsPerGene:
    def test_basic_collection(self, config, sample_group):
        result = collect_locusids_per_gene(sample_group, config)
        assert result["coi"] == "LOC1"
        assert result["16s"] == "LOC2"
        assert result["18s"] == "LOC3"

    def test_empty_group(self, config):
        empty = pd.DataFrame({
            "species_voucher_new": [],
            "organism": [],
            "gene_type": [],
            "LocusID": [],
        })
        result = collect_locusids_per_gene(empty, config)
        for gene in config.gene_order:
            assert result[gene] == ""

    def test_multi_gene_row(self, config):
        group = pd.DataFrame({
            "gene_type": ["coi,16s"],
            "LocusID": ["LOC1"],
        })
        result = collect_locusids_per_gene(group, config)
        assert "LOC1" in result["coi"]
        assert "LOC1" in result["16s"]

    def test_deduplication(self, config):
        group = pd.DataFrame({
            "gene_type": ["coi", "coi"],
            "LocusID": ["LOC1", "LOC1"],
        })
        result = collect_locusids_per_gene(group, config)
        assert result["coi"] == "LOC1"

    def test_nan_gene_type_skipped(self, config):
        group = pd.DataFrame({
            "gene_type": [np.nan],
            "LocusID": ["LOC1"],
        })
        result = collect_locusids_per_gene(group, config)
        assert result["coi"] == ""

    def test_mtgenome_propagation(self, config):
        group = pd.DataFrame({
            "gene_type": ["mtgenome"],
            "LocusID": ["MG1"],
        })
        result = collect_locusids_per_gene(group, config)
        assert "MG1" in result["coi"]
        assert "MG1" in result["16s"]


# --- metadata functions ---


class TestMetadataFirstNonempty:
    def test_basic(self, sample_group):
        result = metadata_first_nonempty(sample_group, ["Class", "Order"])
        assert result["Class"] == "Mammalia"
        assert result["Order"] == "Primates"

    def test_missing_column(self, sample_group):
        result = metadata_first_nonempty(sample_group, ["Class", "nonexistent"])
        assert result["nonexistent"] == ""

    def test_first_nonempty_skips_nan(self):
        df = pd.DataFrame({"col": [np.nan, "val"]})
        result = metadata_first_nonempty(df, ["col"])
        assert result["col"] == "val"


class TestMetadataAll:
    def test_unique_values(self):
        df = pd.DataFrame({"col": ["a", "b", "c"]})
        result = metadata_all(df, ["col"])
        assert result["col"] == "a;b;c"

    def test_deduplicates(self):
        df = pd.DataFrame({"col": ["a", "a", "b"]})
        result = metadata_all(df, ["col"])
        assert result["col"] == "a;b"

    def test_missing_column(self):
        df = pd.DataFrame({"a": [1]})
        result = metadata_all(df, ["nonexistent"])
        assert result["nonexistent"] == ""

    def test_skips_nan_and_empty(self):
        df = pd.DataFrame({"col": ["a", np.nan, "", "b"]})
        result = metadata_all(df, ["col"])
        assert result["col"] == "a;b"

    def test_semicolon_split_values(self):
        df = pd.DataFrame({"col": ["a;b", "c"]})
        result = metadata_all(df, ["col"])
        assert "a" in result["col"]
        assert "b" in result["col"]
        assert "c" in result["col"]


class TestMetadataPerGene:
    def test_basic(self):
        df = pd.DataFrame({
            "gene_type": ["coi", "16s"],
            "country": ["China", "Japan"],
        })
        result = metadata_per_gene(df, ["coi", "16s"], ["country"])
        assert "China" in result["country"]
        assert "Japan" in result["country"]

    def test_missing_column(self):
        df = pd.DataFrame({"gene_type": ["coi"], "nonexistent": ["x"]})
        result = metadata_per_gene(df, ["coi"], ["missing"])
        assert result["missing"] == ""


class TestMetadataAdvanced:
    def test_basic(self):
        df = pd.DataFrame({
            "LocusID": ["LOC1", "LOC2"],
            "gene_type": ["coi", "16s"],
            "country": ["China", "Japan"],
        })
        gene_locusids = {"coi": "LOC1", "16s": "LOC2"}
        result = metadata_advanced(df, gene_locusids, ["country"])
        assert result["country"] != ""

    def test_empty_data(self):
        df = pd.DataFrame({
            "LocusID": [],
            "gene_type": [],
            "country": [],
        })
        result = metadata_advanced(df, {"coi": ""}, ["country"])
        assert result["country"] == ""


# --- process_group ---


class TestProcessGroup:
    def test_basic_group(self, config, sample_group):
        result = process_group(sample_group, config)
        assert result["species_voucher_new"] == "sp1"
        assert result["organism"] == "Homo sapiens"
        assert result["coi"] == "LOC1"
        assert result["16s"] == "LOC2"
        assert result["18s"] == "LOC3"

    def test_metadata_mode_first_nonempty(self, sample_group):
        config = OrganizeConfig(metadata_mode="first_nonempty")
        result = process_group(sample_group, config)
        assert result["Class"] == "Mammalia"


# --- build_header ---


class TestBuildHeader:
    def test_basic_order(self, config):
        sample_row = {"species_voucher_new": "", "organism": "",
                      "coi": "", "16s": "", "18s": ""}
        header = build_header(config, sample_row)
        assert header[0] == "species_voucher_new"
        assert header[1] == "organism"
        # gene_order comes next
        assert "coi" in header
        assert "16s" in header

    def test_no_duplicates(self, config):
        sample_row = {"species_voucher_new": "", "organism": "",
                      "coi": "", "TaxonID": ""}
        header = build_header(config, sample_row)
        assert len(header) == len(set(header))


# --- organize() integration ---


class TestOrganizeIntegration:
    def test_basic_organize(self, tmp_path):
        input_csv = tmp_path / "input.csv"
        input_csv.write_text(
            "species_voucher_new,organism,gene_type,LocusID,Class\n"
            "sp1,Homo sapiens,coi,LOC1,Mammalia\n"
            "sp1,Homo sapiens,16s,LOC2,Mammalia\n"
            "sp2,Mus musculus,coi,LOC3,Mammalia\n"
        )
        output_csv = tmp_path / "output.csv"
        config = OrganizeConfig()
        result = organize(str(input_csv), str(output_csv), config)
        assert result.success is True
        assert result.rows == 2
        assert output_csv.exists()

        df = pd.read_csv(output_csv)
        assert "species_voucher_new" in df.columns
        assert "coi" in df.columns
        assert len(df) == 2

    def test_organize_missing_mandatory_column(self, tmp_path):
        input_csv = tmp_path / "input.csv"
        input_csv.write_text("a,b\n1,2\n")
        output_csv = tmp_path / "output.csv"
        result = organize(str(input_csv), str(output_csv))
        assert result.success is False

    def test_organize_creates_output_dir(self, tmp_path):
        input_csv = tmp_path / "input.csv"
        input_csv.write_text(
            "species_voucher_new,organism,gene_type,LocusID\n"
            "sp1,Org1,coi,LOC1\n"
        )
        output_path = tmp_path / "subdir" / "output.csv"
        result = organize(str(input_csv), str(output_path))
        assert result.success is True
        assert output_path.exists()

    def test_organize_with_nan_group(self, tmp_path):
        input_csv = tmp_path / "input.csv"
        input_csv.write_text(
            "species_voucher_new,organism,gene_type,LocusID\n"
            ",Org1,coi,LOC1\n"
            "sp1,Org2,16s,LOC2\n"
        )
        output_csv = tmp_path / "output.csv"
        result = organize(str(input_csv), str(output_csv))
        assert result.success is True
        df = pd.read_csv(output_csv)
        assert len(df) == 2
