import argparse
import os

import numpy as np
import pandas as pd
import pytest

from g2t.utils import (
    StepResult,
    build_col_lookup,
    col_key,
    is_valid_output,
    length_in_range,
    parse_interval,
    read_table,
    resolve_col,
    safe_concat,
    sanitize_string,
    str2bool,
    write_csv,
)


# --- StepResult ---


class TestStepResult:
    def test_defaults(self):
        r = StepResult(success=True, output_file="out.csv")
        assert r.rows == 0
        assert r.elapsed == 0.0

    def test_full(self):
        r = StepResult(success=False, output_file="err.csv", rows=42, elapsed=1.5)
        assert not r.success
        assert r.rows == 42


# --- str2bool ---


class TestStr2Bool:
    @pytest.mark.parametrize("val", ["yes", "true", "t", "1", "y", "on", "YES", "True", "ON"])
    def test_truthy(self, val):
        assert str2bool(val) is True

    @pytest.mark.parametrize("val", ["no", "false", "f", "0", "n", "off", "NO", "False", "OFF"])
    def test_falsy(self, val):
        assert str2bool(val) is False

    def test_bool_passthrough(self):
        assert str2bool(True) is True
        assert str2bool(False) is False

    def test_invalid_raises(self):
        with pytest.raises(argparse.ArgumentTypeError):
            str2bool("maybe")


# --- sanitize_string ---


class TestSanitizeString:
    def test_none(self):
        assert sanitize_string(None) == ""

    def test_nan(self):
        assert sanitize_string(float("nan")) == ""

    def test_empty(self):
        assert sanitize_string("") == ""
        assert sanitize_string("   ") == ""

    def test_basic(self):
        assert sanitize_string("Hello World") == "Hello_World"

    def test_quotes_removed(self):
        assert sanitize_string("it's \"quoted\"") == "its_quoted"

    def test_parens_removed(self):
        assert sanitize_string("name (copy)") == "name_copy"

    def test_fullwidth_space(self):
        assert sanitize_string("name　value") == "name_value"

    def test_special_chars_replaced(self):
        result = sanitize_string("a<b>c:d/e")
        assert "<" not in result
        assert ">" not in result
        assert ":" not in result
        assert "/" not in result

    def test_underscore_consolidation(self):
        assert sanitize_string("a___b") == "a_b"

    def test_strip_trailing_underscores(self):
        assert sanitize_string("_hello_") == "hello"


# --- safe_concat ---


class TestSafeConcat:
    def test_both_nonempty(self):
        assert safe_concat("a", "b") == "a_b"

    def test_first_empty(self):
        assert safe_concat("", "b") == "b"

    def test_second_empty(self):
        assert safe_concat("a", "") == "a"

    def test_both_empty(self):
        assert safe_concat("", "") == ""

    def test_none_handling(self):
        assert safe_concat(None, "b") == "b"
        assert safe_concat("a", None) == "a"


# --- read_table ---


class TestReadTable:
    def test_csv(self, tmp_path):
        f = tmp_path / "test.csv"
        pd.DataFrame({"a": ["1"], "b": ["2"]}).to_csv(f, index=False)
        df = read_table(str(f))
        assert list(df.columns) == ["a", "b"]
        assert len(df) == 1

    def test_tsv(self, tmp_path):
        f = tmp_path / "test.tsv"
        pd.DataFrame({"a": ["1"], "b": ["2"]}).to_csv(f, index=False, sep="\t")
        df = read_table(str(f))
        assert len(df) == 1

    def test_txt_as_tsv(self, tmp_path):
        f = tmp_path / "test.txt"
        pd.DataFrame({"x": ["val"]}).to_csv(f, index=False, sep="\t")
        df = read_table(str(f))
        assert len(df) == 1


# --- write_csv ---


class TestWriteCsv:
    def test_creates_file(self, tmp_path):
        out = tmp_path / "sub" / "out.csv"
        df = pd.DataFrame({"col": ["val"]})
        write_csv(df, str(out))
        assert out.exists()
        loaded = pd.read_csv(str(out))
        assert len(loaded) == 1

    def test_creates_parent_dir(self, tmp_path):
        out = tmp_path / "deep" / "nested" / "out.csv"
        write_csv(pd.DataFrame({"a": [1]}), str(out))
        assert out.exists()


# --- parse_interval ---


class TestParseInterval:
    def test_colon_range(self):
        assert parse_interval("100:500") == (100, 500)

    def test_comma_range(self):
        assert parse_interval("100,500") == (100, 500)

    def test_single_value(self):
        lo, hi = parse_interval("100")
        assert lo == 100
        assert hi is None

    def test_empty(self):
        assert parse_interval("") == (None, None)

    def test_partial_low(self):
        lo, hi = parse_interval("100:")
        assert lo == 100
        assert hi is None

    def test_partial_high(self):
        lo, hi = parse_interval(":500")
        assert lo is None
        assert hi == 500


# --- length_in_range ---


class TestLengthInRange:
    def test_within(self):
        assert length_in_range(300, "100:500") is True

    def test_below(self):
        assert length_in_range(50, "100:500") is False

    def test_above(self):
        assert length_in_range(600, "100:500") is False

    def test_none_string(self):
        assert length_in_range(999, "none") is True

    def test_empty(self):
        assert length_in_range(999, "") is True

    def test_non_numeric_length(self):
        assert length_in_range("abc", "100:500") is True


# --- is_valid_output ---


class TestIsValidOutput:
    def test_existing_nonempty(self, tmp_path):
        f = tmp_path / "file.csv"
        f.write_text("data")
        assert is_valid_output(str(f)) is True

    def test_missing(self, tmp_path):
        assert is_valid_output(str(tmp_path / "nope.csv")) is False

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.csv"
        f.write_text("")
        assert is_valid_output(str(f)) is False

    def test_directory(self, tmp_path):
        assert is_valid_output(str(tmp_path)) is False


# --- col_key ---


class TestColKey:
    def test_lowercase(self):
        assert col_key("Hello") == "hello"

    def test_spaces_to_underscores(self):
        assert col_key("Gene Type") == "gene_type"

    def test_already_normalized(self):
        assert col_key("gene_type") == "gene_type"


# --- resolve_col ---


class TestResolveCol:
    def test_exact_match(self):
        df = pd.DataFrame(columns=["Gene", "Length"])
        assert resolve_col(df, "Gene") == "Gene"

    def test_case_insensitive(self):
        df = pd.DataFrame(columns=["gene", "Length"])
        assert resolve_col(df, "GENE") == "gene"

    def test_space_to_underscore(self):
        df = pd.DataFrame(columns=["Gene Type", "Length"])
        assert resolve_col(df, "gene_type") == "Gene Type"

    def test_underscore_to_space(self):
        df = pd.DataFrame(columns=["Gene Type", "Length"])
        assert resolve_col(df, "Gene_Type") == "Gene Type"

    def test_missing_returns_none(self):
        df = pd.DataFrame(columns=["Gene", "Length"])
        assert resolve_col(df, "Missing") is None

    def test_with_lookup(self):
        df = pd.DataFrame(columns=["Gene", "Length"])
        lookup = build_col_lookup(df)
        assert resolve_col(df, "gene", lookup) == "Gene"


# --- build_col_lookup ---


class TestBuildColLookup:
    def test_basic(self):
        df = pd.DataFrame(columns=["Gene Type", "Length"])
        lookup = build_col_lookup(df)
        assert lookup["gene_type"] == "Gene Type"
        assert lookup["length"] == "Length"
