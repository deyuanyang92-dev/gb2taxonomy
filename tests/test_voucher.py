"""Tests for voucher.py functions."""

import os
import re

import pandas as pd
import pytest

from g2t.voucher import (
    build_species_vouchers,
    log,
    normalize_columns_inplace,
    now_ts,
)


class TestNowTs:
    def test_format(self):
        ts = now_ts()
        assert re.match(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", ts)


class TestLog:
    def test_prints_to_stdout(self, capsys):
        log("test message", quiet=False)
        captured = capsys.readouterr()
        assert "test message" in captured.out

    def test_quiet_suppresses_stdout(self, capsys):
        log("hidden", quiet=True)
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_writes_to_file(self, tmp_path):
        log_file = tmp_path / "test.log"
        log("file message", quiet=True, log_file_path=str(log_file))
        assert log_file.exists()
        assert "file message" in log_file.read_text()

    def test_no_file_when_none(self, tmp_path):
        log("no file", quiet=True, log_file_path=None)
        # Should not crash


class TestNormalizeColumnsInplace:
    def test_basic_normalization(self):
        df = pd.DataFrame({"First Name": [1], "Last Name": [2]})
        normalize_columns_inplace(df)
        assert "first_name" in df.columns
        assert "last_name" in df.columns

    def test_deduplication(self):
        df = pd.DataFrame({"col": [1], "col": [2]})
        # When two columns have the same name after normalization
        df = pd.DataFrame({"Col": [1], "col": [2]})
        normalize_columns_inplace(df)
        assert df.columns[0] == "col"
        assert "dup" in df.columns[1]

    def test_lowercases_and_replaces_spaces(self):
        df = pd.DataFrame({"First Name": [1]})
        normalize_columns_inplace(df)
        assert "first_name" in df.columns

    def test_preserves_data(self):
        df = pd.DataFrame({"A B": [42]})
        normalize_columns_inplace(df)
        assert df.iloc[0, 0] == 42


class TestBuildSpeciesVouchersIntegration:
    @pytest.fixture
    def input_csv(self, tmp_path):
        """Create a minimal input CSV for voucher testing."""
        csv_path = tmp_path / "input.csv"
        csv_path.write_text(
            "organism,gene_type,LocusID,ACCESSION,specimen_voucher,isolate\n"
            "Homo sapiens,coi,LOC1,ACC1,V001,\n"
            "Homo sapiens,16s,LOC2,ACC2,V001,\n"
            "Mus musculus,coi,LOC3,ACC3,,ISO1\n"
            "Canis lupus,18s,LOC4,ACC4,,\n"
        )
        return str(csv_path)

    def test_basic_build(self, tmp_path, input_csv):
        out_dir = tmp_path / "output"
        result = build_species_vouchers(
            file_path=input_csv,
            output_dir=str(out_dir),
            quiet=True,
        )
        assert result.success is True
        assert result.rows == 4
        assert os.path.exists(os.path.join(str(out_dir), "updated_species_voucher.csv"))

    def test_output_has_species_voucher_new(self, tmp_path, input_csv):
        out_dir = tmp_path / "output"
        build_species_vouchers(
            file_path=input_csv,
            output_dir=str(out_dir),
            quiet=True,
        )
        df = pd.read_csv(os.path.join(str(out_dir), "updated_species_voucher.csv"))
        assert "species_voucher_new" in df.columns

    def test_count_file_created(self, tmp_path, input_csv):
        out_dir = tmp_path / "output"
        build_species_vouchers(
            file_path=input_csv,
            output_dir=str(out_dir),
            quiet=True,
        )
        count_file = os.path.join(str(out_dir), "count_species_voucher.txt")
        assert os.path.exists(count_file)

    def test_pse_only_mode(self, tmp_path, input_csv):
        out_dir = tmp_path / "output"
        result = build_species_vouchers(
            file_path=input_csv,
            output_dir=str(out_dir),
            combine_order="pse_only",
            quiet=True,
        )
        assert result.success is True
        df = pd.read_csv(os.path.join(str(out_dir), "updated_species_voucher.csv"))
        assert "species_voucher_new" in df.columns

    def test_simple_csv_generated(self, tmp_path, input_csv):
        out_dir = tmp_path / "output"
        build_species_vouchers(
            file_path=input_csv,
            output_dir=str(out_dir),
            if_generate_simple_csv=True,
            quiet=True,
        )
        simple_path = os.path.join(str(out_dir), "shorten_species_voucher.csv")
        assert os.path.exists(simple_path)

    def test_normalize_columns(self, tmp_path):
        csv_path = tmp_path / "input.csv"
        csv_path.write_text(
            "organism,gene type,Locus ID,ACCESSION\n"
            "Org1,coi,L1,A1\n"
        )
        out_dir = tmp_path / "output"
        result = build_species_vouchers(
            file_path=str(csv_path),
            output_dir=str(out_dir),
            normalize_column_names=True,
            quiet=True,
        )
        assert result.success is True

    def test_fallback_to_accession(self, tmp_path):
        csv_path = tmp_path / "input.csv"
        csv_path.write_text(
            "organism,gene_type,LocusID,ACCESSION\n"
            "Org1,coi,L1,ACC_FALLBACK\n"
        )
        out_dir = tmp_path / "output"
        build_species_vouchers(
            file_path=str(csv_path),
            output_dir=str(out_dir),
            quiet=True,
        )
        df = pd.read_csv(os.path.join(str(out_dir), "updated_species_voucher.csv"))
        # ACCESSION is the fallback, so species_voucher_pse should have ACC_FALLBACK
        assert "ACC_FALLBACK" in str(df["species_voucher_pse"].values)
