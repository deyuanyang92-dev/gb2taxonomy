"""Voucher standardisation and metadata curation (g2t.curate)."""

import pandas as pd
import pytest

from g2t.curate import apply_updates, make_template, normalize_voucher, standardize_vouchers


@pytest.mark.parametrize("raw,norm", [
    ("COI_ZMMU_MSU_WS12387_XZ5022", "ZMMU_MSU_WS12387_XZ5022"),
    ("28S_ZMMU_MSU_WS12387_XZ5022", "ZMMU_MSU_WS12387_XZ5022"),
    ("ZMMU:WS30980", "ZMMU_WS30980"),
    ("ZMMU MSU WS14906", "ZMMU_MSU_WS14906"),
    ("28S_ZMMU_WS3020", "ZMMU_WS3020"),
    ("WS0397", "WS0397"),
    ("07PROBE-05367", "07PROBE-05367"),
    ("BNSB0286", "BNSB0286"),
    ("4", "4"),
    ("", ""),
])
def test_normalize_voucher(raw, norm):
    assert normalize_voucher(raw) == norm


@pytest.mark.parametrize("variants,expected,consistent", [
    (["ZMMU:WS30980", "ZMMU_WS30980"], "ZMMU_WS30980", True),
    (["COI_ZMMU_MSU_WS12387_XZ5022", "28S_ZMMU_MSU_WS12387_XZ5022"], "ZMMU_MSU_WS12387_XZ5022", True),
    (["COI_ZMMU_MSU_WS2585", "28S_ZMMU_WS2585", "WS2585"], "ZMMU_MSU_WS2585", True),
    (["28S_ZMMU_MSU_WS16798_Pr20", "Pr20"], "ZMMU_MSU_WS16798_Pr20", True),
    (["COI_ZMMU_MSU_WS3017", "WS3017"], "ZMMU_MSU_WS3017", True),
    (["WS399", "WS400"], "WS399", False),
    ([], "", True),
    (["", "  "], "", True),
])
def test_standardize_vouchers(variants, expected, consistent):
    canon, ok = standardize_vouchers(variants)
    assert canon == expected and ok is consistent


@pytest.fixture
def matrix():
    return pd.DataFrame([
        {"species_voucher_new": "Pc_WS399", "organism": "Priapulus caudatus", "coi": "ON792938.1",
         "28s": "ON793000.1", "16s": "OP279987.1", "voucher_standardized": "ZMMU_MSU_WS399",
         "lat_lon": "", "geo_loc_name": "", "Ref1Title": "Direct Submission"},
        {"species_voucher_new": "Pc_WS3020", "organism": "Priapulus caudatus", "coi": "ON792926.1",
         "28s": "", "16s": "OP279978.1", "voucher_standardized": "ZMMU_MSU_WS3020",
         "lat_lon": "66.55 N 33.10 E", "geo_loc_name": "Russia", "Ref1Title": ""},
        {"species_voucher_new": "Hs_WS3020", "organism": "Halicryptus spinulosus", "coi": "",
         "28s": "ON792990.1", "16s": "", "voucher_standardized": "ZMMU_WS3020",
         "lat_lon": "", "geo_loc_name": "", "Ref1Title": ""},
    ])


GENES = ["coi", "28s", "16s"]


def test_update_by_accession_with_or_without_version(matrix):
    upd = pd.DataFrame([{"accession": "ON792938", "lat_lon": "66.55 N 33.10 E"},
                        {"accession": "OP279978.1", "geo_loc_name": "Russia: White Sea"}])
    out, log, problems = apply_updates(matrix, upd, gene_columns=GENES)
    assert out.loc[0, "lat_lon"] == "66.55 N 33.10 E"
    assert out.loc[1, "geo_loc_name"] == "Russia: White Sea"
    assert set(log["field"]) == {"lat_lon", "geo_loc_name"}
    assert log.loc[log.field == "geo_loc_name", "old_value"].iloc[0] == "Russia"
    assert problems.empty
    assert matrix.loc[1, "geo_loc_name"] == "Russia"            # input not modified


def test_update_by_voucher_any_written_form(matrix):
    upd = pd.DataFrame([{"voucher": "COI_ZMMU_MSU_WS399", "organism": "Priapulus tuberculatospinosus",
                         "Ref1Title": "Population genetics reveals ..."}])
    out, log, problems = apply_updates(matrix, upd, gene_columns=GENES)
    assert out.loc[0, "organism"] == "Priapulus tuberculatospinosus"
    assert out.loc[0, "curated_fields"] == "organism; Ref1Title"
    assert len(log) == 2 and problems.empty


def test_empty_cells_do_not_overwrite(matrix):
    upd = pd.DataFrame([{"accession": "ON792926", "lat_lon": "", "geo_loc_name": "Norway"}])
    out, log, _ = apply_updates(matrix, upd, gene_columns=GENES)
    assert out.loc[1, "lat_lon"] == "66.55 N 33.10 E" and out.loc[1, "geo_loc_name"] == "Norway"
    assert list(log["field"]) == ["geo_loc_name"]


def test_no_match_and_unknown_column_reported(matrix):
    upd = pd.DataFrame([{"accession": "XX000001", "lat_lon": "1 N 1 E"},
                        {"accession": "ON792938", "new_field": "x"}])
    out, log, problems = apply_updates(matrix, upd, gene_columns=GENES)
    assert "no matrix row" in " ".join(problems["problem"])
    assert out.loc[0, "new_field"] == "x"                      # unknown columns are added
    assert "new column" in " ".join(problems["problem"])


def test_ambiguous_voucher_is_not_applied(matrix):
    upd = pd.DataFrame([{"voucher": "WS3020", "lat_lon": "1 N 1 E"}])   # Pc_WS3020 and Hs_WS3020
    out, log, problems = apply_updates(matrix, upd, gene_columns=GENES)
    assert log.empty and "2 matrix rows" in problems["problem"].iloc[0]
    upd2 = pd.DataFrame([{"voucher": "WS3020", "organism_match": "Halicryptus spinulosus", "lat_lon": "1 N 1 E"}])
    out2, log2, problems2 = apply_updates(matrix, upd2, gene_columns=GENES)
    assert out2.loc[2, "lat_lon"] == "1 N 1 E" and problems2.empty


def test_accession_and_voucher_must_agree(matrix):
    upd = pd.DataFrame([{"accession": "ON792938", "voucher": "WS3020", "lat_lon": "1 N 1 E"}])
    _, log, problems = apply_updates(matrix, upd, gene_columns=GENES)
    assert log.empty and "point to different rows" in problems["problem"].iloc[0]


def test_template_lists_specimens_and_fields(matrix):
    t = make_template(matrix, gene_columns=GENES)
    assert list(t.columns[:3]) == ["voucher", "accession", "organism_match"]
    assert len(t) == 3
    assert t.loc[0, "accession"] == "ON792938.1"
    assert "lat_lon" in t.columns and "organism" in t.columns
