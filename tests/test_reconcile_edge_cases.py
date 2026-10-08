"""Edge cases found in review (R1-R8): key collisions, transitive merges, antimeridian, cores."""

import pandas as pd
import pytest

from g2t.reconcile import reconcile_dataframe, voucher_cores


def row(acc, key, voucher, gene, org="Priapulus caudatus", **kw):
    r = dict(ACCESSION=acc, LocusID=acc, organism=org, gene_type=gene, specimen_voucher=voucher,
             species_voucher_new=key, collection_date="", lat_lon="", geo_loc_name="", collected_by="",
             Ref1Title="", Ref1Authors="")
    r.update(kw)
    return r


def specimens(out):
    return out.groupby("species_voucher_new")["ACCESSION"].apply(lambda s: tuple(sorted(s))).tolist()


def test_r1_three_specimens_with_one_core_stay_three():
    rows = []
    for yr, num in ((2012, "77"), (2013, "88"), (2014, "99")):
        d = f"{yr}-06-01"
        rows.append(row(f"C{yr}", f"Pc_COI_{yr}", f"WS12{num}", "coi", collection_date=d))
        rows.append(row(f"N{yr}", f"Pc_28S_{yr}", f"WS12{num}", "28s", collection_date=d))
    # force a shared core for all three by using the same voucher with different dates
    for r in rows:
        r["specimen_voucher"] = "WS1234"
    out, _ = reconcile_dataframe(pd.DataFrame(rows))
    assert sorted(specimens(out)) == [("C2012", "N2012"), ("C2013", "N2013"), ("C2014", "N2014")]
    assert out["species_voucher_new"].nunique() == 3


def test_r2_no_transitive_bypass_of_same_gene_rule():
    a = row("A1", "kA", "WS399", "coi", collection_date="2019-06-12", collected_by="Smith",
            geo_loc_name="Russia: White Sea")
    b = row("B1", "kB", "WS399", "28s", collection_date="2019-06-12", lat_lon="66.5 N 33.1 E")
    c = row("C1", "kC", "WS399", "coi", lat_lon="66.5 N 33.1 E", collected_by="Smith",
            geo_loc_name="Russia: White Sea")
    out, rep = reconcile_dataframe(pd.DataFrame([a, b, c]))
    groups = specimens(out)
    assert len(groups) == 2
    assert not any(set(g) >= {"A1", "C1"} for g in groups)
    assert rep["decision"].str.contains("same gene twice").any()


def test_r3_conflict_across_different_cores_blocks_merge():
    a = row("A1", "kA", "WS399", "coi", collection_date="2001-01-01", geo_loc_name="Norway: Oslo", Ref1Title="Paper X")
    b1 = row("B1", "kB", "WS399", "28s", Ref1Title="Paper X")
    b2 = row("B2", "kB", "WS500", "18s", Ref1Title="Paper X")
    c = row("C1", "kC", "WS500", "h3", collection_date="2019-06-12", geo_loc_name="Russia: White Sea",
            Ref1Title="Paper X")
    out, rep = reconcile_dataframe(pd.DataFrame([a, b1, b2, c]))
    assert not any(set(g) >= {"A1", "C1"} for g in specimens(out))
    assert rep["decision"].str.contains("conflicting groups").any()


def test_r4_antimeridian_is_close():
    a = row("A1", "kA", "WS399", "coi", lat_lon="50.0 N 179.99 E", collected_by="Smith",
            geo_loc_name="Russia: Kuril")
    b = row("B1", "kB", "WS399", "28s", lat_lon="50.0 N 179.99 W", collected_by="Smith",
            geo_loc_name="Russia: Kuril")
    out, rep = reconcile_dataframe(pd.DataFrame([a, b]))
    assert out["species_voucher_new"].nunique() == 1
    assert not rep["evidence"].str.contains("conflict").any()


@pytest.mark.parametrize("raw,expected", [
    ("MNHN-IU-2014-1234", ["IU20141234", "#20141234"]),
    ("WS399 (holotype)", ["WS399"]),
    ("WS399.", ["WS399"]),
    ("WS399, WS400", ["WS399", "WS400"]),
    ("WS399/WS400", ["WS399", "WS400"]),
    ("USNM 123456", ["USNM123456", "#123456"]),
    ("123456", ["123456", "#123456"]),
    ("2019-123", ["2019123", "#2019123"]),
])
def test_r5_voucher_cores(raw, expected):
    assert voucher_cores(raw) == expected


def test_r5_years_separate_specimens_and_bare_numbers_meet():
    a = row("A1", "kA", "MNHN-IU-2013-1234", "coi", Ref1Title="Paper X")
    b = row("B1", "kB", "MNHN-IU-2014-1234", "28s", Ref1Title="Paper X")
    c = row("C1", "kC", "USNM 123456", "coi", Ref1Title="Paper Y")
    d = row("D1", "kD", "123456", "28s", Ref1Title="Paper Y")
    out, _ = reconcile_dataframe(pd.DataFrame([a, b, c, d]))
    groups = specimens(out)
    assert ("A1",) in groups and ("B1",) in groups
    assert ("C1", "D1") in groups


def test_r6_missing_organism_column_is_an_error():
    a = row("A1", "kA", "WS399", "coi", Ref1Title="Paper X")
    b = row("B1", "kB", "WS399", "28s", org="Halicryptus spinulosus", Ref1Title="Paper X")
    with pytest.raises(ValueError, match="organism"):
        reconcile_dataframe(pd.DataFrame([a, b]).drop(columns=["organism"]))


def test_r6_empty_organism_not_merged():
    a = row("A1", "kA", "WS399", "coi", org="", Ref1Title="Paper X")
    b = row("B1", "kB", "WS399", "28s", org="", Ref1Title="Paper X")
    out, rep = reconcile_dataframe(pd.DataFrame([a, b]))
    assert out["species_voucher_new"].nunique() == 2
    assert (rep["decision"] == "not merged: organism missing").all()


def test_r7_multi_gene_types_count_as_overlap():
    a = row("A1", "kA", "WS399", "coi,16s", collected_by="Smith", geo_loc_name="Russia: White Sea")
    b = row("B1", "kB", "WS399", "16s", collected_by="Smith", geo_loc_name="Russia: White Sea")
    out, rep = reconcile_dataframe(pd.DataFrame([a, b]))
    assert out["species_voucher_new"].nunique() == 2
    assert rep["decision"].str.contains("same gene").any()


def test_r8_clone_before_strain():
    a = row("A1", "kA", "", "coi", clone="C1234", strain="S5678", Ref1Title="P")
    out, _ = reconcile_dataframe(pd.DataFrame([a]))
    assert out.loc[0, "voucher_core"] == "C1234"
