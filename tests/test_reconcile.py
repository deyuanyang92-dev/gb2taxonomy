"""Tests for reconcile.py: evidence-based merging of voucher variants."""

import pandas as pd
import pytest

from g2t.reconcile import (
    ReconcileConfig,
    compare_groups,
    reconcile_dataframe,
    voucher_core,
)

PAPER = "Population genetics reveals a complex Priapulus caudatus"


def rec(acc, org, gene, voucher, key, **kw):
    row = dict(ACCESSION=acc, organism=org, gene_type=gene, specimen_voucher=voucher, isolate="",
               species_voucher_new=key, Ref1Title="", Ref2Title="", Ref1Authors="", collection_date="",
               lat_lon="", geo_loc_name="", country="", collected_by="")
    row.update(kw)
    return row


class TestVoucherCore:
    @pytest.mark.parametrize("raw,core", [
        ("COI_ZMMU_MSU_WS399", "WS399"),
        ("28S_ZMMU_WS399", "WS399"),
        ("WS399", "WS399"),
        ("ZMMU WS-399", "WS399"),
        ("WS399_COI", "WS399"),
        ("18S-XZ4579", "XZ4579"),
        ("ZMH V13492", "V13492"),
        ("", ""),
        ("A1", ""),            # too short to be a reliable identifier
    ])
    def test_core(self, raw, core):
        assert voucher_core(raw) == core


class TestCompare:
    def test_shared_paper_is_strong(self):
        a = pd.DataFrame([rec("A", "P c", "coi", "COI_WS399", "k1", Ref1Title=PAPER)])
        b = pd.DataFrame([rec("B", "P c", "28s", "28S_WS399", "k2", Ref1Title=PAPER)])
        ev = compare_groups(a, b)
        assert ev.strong and not ev.conflicts

    def test_direct_submission_not_evidence(self):
        a = pd.DataFrame([rec("A", "P c", "coi", "WS399", "k1", Ref1Title="Direct Submission")])
        b = pd.DataFrame([rec("B", "P c", "28s", "WS399x", "k2", Ref1Title="Direct Submission")])
        ev = compare_groups(a, b)
        assert not ev.strong and not ev.moderate

    def test_date_conflict(self):
        a = pd.DataFrame([rec("A", "P c", "coi", "WS399", "k1", Ref1Title=PAPER, collection_date="12-Jun-2019")])
        b = pd.DataFrame([rec("B", "P c", "28s", "WS399", "k2", Ref1Title=PAPER, collection_date="03-Aug-2015")])
        assert compare_groups(a, b).conflicts

    def test_partial_date_compatible(self):
        a = pd.DataFrame([rec("A", "P c", "coi", "WS399", "k1", collection_date="2019")])
        b = pd.DataFrame([rec("B", "P c", "28s", "WS399", "k2", collection_date="12-Jun-2019")])
        ev = compare_groups(a, b)
        assert not ev.conflicts

    def test_country_conflict(self):
        a = pd.DataFrame([rec("A", "P c", "coi", "WS399", "k1", geo_loc_name="Russia: White Sea")])
        b = pd.DataFrame([rec("B", "P c", "28s", "WS399", "k2", geo_loc_name="Norway")])
        assert compare_groups(a, b).conflicts

    def test_latlon_conflict_and_match(self):
        a = pd.DataFrame([rec("A", "P c", "coi", "WS399", "k1", lat_lon="66.55 N 33.10 E")])
        b = pd.DataFrame([rec("B", "P c", "28s", "WS399", "k2", lat_lon="66.55 N 33.10 E")])
        c = pd.DataFrame([rec("C", "P c", "18s", "WS399", "k3", lat_lon="10.00 S 33.10 E")])
        assert compare_groups(a, b).strong
        assert compare_groups(a, c).conflicts


class TestReconcile:
    def test_merges_gene_prefixed_vouchers_with_evidence(self):
        df = pd.DataFrame([
            rec("A", "Priapulus caudatus", "coi", "COI_ZMMU_MSU_WS399", "Pc_COI_ZMMU_MSU_WS399", Ref1Title=PAPER),
            rec("B", "Priapulus caudatus", "28s", "28S_ZMMU_WS399", "Pc_28S_ZMMU_WS399", Ref1Title=PAPER),
            rec("C", "Priapulus caudatus", "16s", "WS399", "Pc_WS399", Ref1Title=PAPER, lat_lon="66.55 N 33.10 E"),
        ])
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1
        assert set(out["match_confidence"]) == {"high"}
        assert set(out["match_basis"]) == {"reconciled"}
        assert (out["species_voucher_g2t"] == df["species_voucher_new"]).all()

    def test_never_merges_across_organisms(self):
        df = pd.DataFrame([
            rec("A", "Halicryptus spinulosus", "28s", "28S_ZMMU_WS3020", "Hs_WS3020", Ref1Title=PAPER),
            rec("B", "Priapulus caudatus", "coi", "COI_ZMMU_MSU_WS3020", "Pc_WS3020", Ref1Title=PAPER),
        ])
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert (report["decision"] == "not merged: different organisms").any()

    def test_no_evidence_not_merged(self):
        df = pd.DataFrame([
            rec("A", "P c", "coi", "COI_WS399", "k1"),
            rec("B", "P c", "28s", "28S_WS399", "k2"),
        ])
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert (report["decision"] == "not merged: insufficient evidence").any()

    def test_conflict_blocks_merge(self):
        df = pd.DataFrame([
            rec("A", "P c", "coi", "COI_WS399", "k1", Ref1Title=PAPER, geo_loc_name="Russia"),
            rec("B", "P c", "28s", "28S_WS399", "k2", Ref1Title=PAPER, geo_loc_name="Chile"),
        ])
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert report["decision"].str.startswith("not merged: conflicting").any()

    def test_same_gene_twice_needs_strong_evidence(self):
        df = pd.DataFrame([
            rec("A", "P c", "coi", "COI_WS399", "k1", collected_by="X", geo_loc_name="Russia: A"),
            rec("B", "P c", "coi", "WS399", "k2", collected_by="X", geo_loc_name="Russia: A"),
        ])
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2

    def test_medium_confidence_threshold(self):
        df = pd.DataFrame([
            rec("A", "P c", "coi", "COI_WS399", "k1", collected_by="X", geo_loc_name="Russia: A"),
            rec("B", "P c", "28s", "WS399", "k2", collected_by="X", geo_loc_name="Russia: A"),
        ])
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1
        assert set(out["match_confidence"]) == {"medium"}
        out2, _ = reconcile_dataframe(df, ReconcileConfig(min_confidence="high"))
        assert out2["species_voucher_new"].nunique() == 2

    def test_untouched_rows_marked_exact(self):
        df = pd.DataFrame([rec("A", "P c", "coi", "V12345", "k1")])
        out, report = reconcile_dataframe(df)
        assert out.loc[0, "match_basis"] == "exact"
        assert out.loc[0, "species_voucher_new"] == "k1"
        assert report.empty
