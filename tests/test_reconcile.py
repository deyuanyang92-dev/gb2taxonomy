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
    row = dict(
        ACCESSION=acc,
        organism=org,
        gene_type=gene,
        specimen_voucher=voucher,
        isolate="",
        species_voucher_new=key,
        Ref1Title="",
        Ref2Title="",
        Ref1Authors="",
        collection_date="",
        lat_lon="",
        geo_loc_name="",
        country="",
        collected_by="",
    )
    row.update(kw)
    return row


class TestVoucherCore:
    @pytest.mark.parametrize(
        "raw,core",
        [
            ("COI_ZMMU_MSU_WS399", "WS399"),
            ("28S_ZMMU_WS399", "WS399"),
            ("WS399", "WS399"),
            ("ZMMU WS-399", "WS399"),
            ("WS399_COI", "WS399"),
            ("18S-XZ4579", "XZ4579"),
            ("ZMH V13492", "V13492"),
            ("", ""),
            ("A1", ""),  # too short to be a reliable identifier
        ],
    )
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
        df = pd.DataFrame(
            [
                rec("A", "Priapulus caudatus", "coi", "COI_ZMMU_MSU_WS399", "Pc_COI_ZMMU_MSU_WS399", Ref1Title=PAPER),
                rec("B", "Priapulus caudatus", "28s", "28S_ZMMU_WS399", "Pc_28S_ZMMU_WS399", Ref1Title=PAPER),
                rec("C", "Priapulus caudatus", "16s", "WS399", "Pc_WS399", Ref1Title=PAPER, lat_lon="66.55 N 33.10 E"),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1
        assert set(out["match_confidence"]) == {"high"}
        assert set(out["match_basis"]) == {"reconciled"}
        assert (out["species_voucher_g2t"] == df["species_voucher_new"]).all()

    def test_never_merges_across_organisms(self):
        df = pd.DataFrame(
            [
                rec("A", "Halicryptus spinulosus", "28s", "28S_ZMMU_WS3020", "Hs_WS3020", Ref1Title=PAPER),
                rec("B", "Priapulus caudatus", "coi", "COI_ZMMU_MSU_WS3020", "Pc_WS3020", Ref1Title=PAPER),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert (report["decision"] == "not merged: different organisms").any()

    def test_no_evidence_not_merged(self):
        df = pd.DataFrame(
            [
                rec("A", "P c", "coi", "COI_WS399", "k1"),
                rec("B", "P c", "28s", "28S_WS399", "k2"),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert (report["decision"] == "not merged: insufficient evidence").any()

    def test_conflict_blocks_merge(self):
        df = pd.DataFrame(
            [
                rec("A", "P c", "coi", "COI_WS399", "k1", Ref1Title=PAPER, geo_loc_name="Russia"),
                rec("B", "P c", "28s", "28S_WS399", "k2", Ref1Title=PAPER, geo_loc_name="Chile"),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert report["decision"].str.startswith("not merged: conflicting").any()

    def test_same_gene_twice_needs_strong_evidence(self):
        df = pd.DataFrame(
            [
                rec("A", "P c", "coi", "COI_WS399", "k1", collected_by="X", geo_loc_name="Russia: A"),
                rec("B", "P c", "coi", "WS399", "k2", collected_by="X", geo_loc_name="Russia: A"),
            ]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2

    def test_medium_confidence_threshold(self):
        df = pd.DataFrame(
            [
                rec("A", "P c", "coi", "COI_WS399", "k1", collected_by="X", geo_loc_name="Russia: A"),
                rec("B", "P c", "28s", "WS399", "k2", collected_by="X", geo_loc_name="Russia: A"),
            ]
        )
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


class TestRowOrderIndependence:
    def test_same_clusters_for_any_row_order(self):
        rows = [
            rec("A", "P c", "coi", "COI_WS399", "kA", Ref1Title=PAPER, Ref2Title="Other paper"),
            rec("B", "P c", "28s", "28S_WS399", "kB", Ref1Title=PAPER, geo_loc_name="Russia"),
            rec("C", "P c", "16s", "16S_WS399", "kC", Ref1Title="Other paper", geo_loc_name="Norway"),
        ]
        results = set()
        for order in ([0, 1, 2], [1, 0, 2], [2, 1, 0], [1, 2, 0]):
            out, _ = reconcile_dataframe(pd.DataFrame([rows[i] for i in order]))
            key = out.set_index("ACCESSION")["species_voucher_new"]
            results.add(frozenset(frozenset(key.index[key == v]) for v in key.unique()))
        assert len(results) == 1

    def test_punctuation_only_voucher_falls_back_to_isolate(self):
        df = pd.DataFrame(
            [
                rec("A", "P c", "coi", "..", "Pc_WS399", isolate="WS399", Ref1Title=PAPER),
                rec("B", "P c", "28s", "28S_ZMMU_WS399", "Pc_28S_ZMMU_WS399", Ref1Title=PAPER),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1


class TestSameVoucherNameVariants:
    """One specimen deposited under slightly different organism names (Polynoidae SIO-BIC A10920)."""

    def test_same_genus_identical_voucher_merged(self):
        df = pd.DataFrame(
            [
                rec("OM007981", "Bathykurila sp.", "18s", "SIO:BIC:A10920", "Bs_A10920"),
                rec("MW783685", "Bathykurila sp.", "mtgenome", "SIO:BIC:A10920", "Bs_A10920"),
                rec("OL960135", "Bathykurila sp. AH-2023", "h3", "SIO:BIC:A10920", "BsAH_A10920"),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1
        assert set(out["match_basis"]) == {"reconciled"}
        assert out["match_evidence"].str.contains("organism names differ").all()

    def test_cf_variant_on_isolate_merged(self):
        df = pd.DataFrame(
            [
                rec("PP968595", "Pelodrilus darlingensis", "12s", "", "Pd_CE17262", isolate="CE17262"),
                rec("PP994587", "Pelodrilus cf. darlingensis", "16s", "", "Pcd_CE17262", isolate="CE17262"),
            ]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1

    def test_identical_code_in_another_field(self):
        """MF001136 has voucher 'personal collection:Christer Erseus:CE14346' and isolate CE14346."""
        df = pd.DataFrame(
            [
                rec(
                    "MF001136",
                    "Capilloventer australis",
                    "coi",
                    "personal collection:Christer Erseus:CE14346",
                    "Ca_pc_CE14346",
                    isolate="CE14346",
                ),
                rec("PP968629", "Capilloventer australis", "12s", "", "Ca_CE14346", isolate="CE14346"),
                rec("PP994526", "Capilloventer australis", "28s", "", "Ca_CE14346", isolate="CE14346"),
            ]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1

    def test_voucher_punctuation_variant_same_organism_merged(self):
        df = pd.DataFrame(
            [
                rec("A", "Eulagisca gigantea", "18s", "SIO:BIC:A6302", "k1"),
                rec("B", "Eulagisca gigantea", "coi", "SIO-BIC A6302", "k2"),
            ]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1

    def test_different_genera_kept_apart_with_warning(self):
        df = pd.DataFrame(
            [
                rec("A", "Halicryptus spinulosus", "28s", "WS3020", "Hs_WS3020"),
                rec("B", "Priapulus caudatus", "coi", "WS3020", "Pc_WS3020"),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert report["decision"].str.startswith("not merged: same voucher, different genera").any()
        assert out["match_evidence"].str.contains("same voucher as").all()

    def test_conflict_still_blocks(self):
        df = pd.DataFrame(
            [
                rec("A", "Bathykurila sp.", "18s", "SIO:BIC:A10920", "k1", geo_loc_name="USA"),
                rec("B", "Bathykurila sp. AH-2023", "h3", "SIO:BIC:A10920", "k2", geo_loc_name="Japan"),
            ]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2


class TestSameGeneSamePaper:
    """Martin et al. 2024 (Haplotaxidae) list CE438 and CE438_2 as two specimens; both have a 28S."""

    def test_shared_paper_alone_does_not_merge_same_gene(self):
        paper = "Towards an integrative revision of Haplotaxidae"
        df = pd.DataFrame(
            [
                rec("PP994571", "Haplotaxidae gen. 1", "28s", "", "k438", isolate="CE438", Ref1Title=paper),
                rec("PP994570", "Haplotaxidae gen. 1", "28s", "", "k438_2", isolate="CE438_2", Ref1Title=paper),
                rec("PP968606", "Haplotaxidae gen. 1", "12s", "", "k438_2", isolate="CE438_2", Ref1Title=paper),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert report["decision"].str.contains("sub-identifier").any()

    def test_same_gene_and_shared_paper_flagged(self):
        paper = "Towards an integrative revision of Haplotaxidae"
        df = pd.DataFrame(
            [
                rec("A", "P c", "28s", "COI_ZMMU_WS399", "k1", Ref1Title=paper),
                rec("B", "P c", "28s", "WS399", "k2", Ref1Title=paper),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert report["decision"].str.contains("same gene in both groups").any()
        assert out["match_evidence"].str.contains("similar voucher").all()

    def test_shared_paper_still_merges_different_genes(self):
        paper = "Towards an integrative revision of Haplotaxidae"
        df = pd.DataFrame(
            [
                rec("A", "P c", "coi", "COI_WS399", "k1", Ref1Title=paper),
                rec("B", "P c", "28s", "28S_WS399", "k2", Ref1Title=paper),
            ]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1


class TestSubIdentifiers:
    """Vouchers that share a core but end in different sub-identifiers are different specimens
    (lot numbers: BIOUG06836-E01/-E02, SIO:BIC:A9392-1/-2, BP2010-044/-045, CE438/CE438_2)."""

    @pytest.mark.parametrize(
        "va,vb",
        [
            ("BIOUG06836-E01", "BIOUG06836-E02"),
            ("SIO:BIC:A9392-1", "SIO:BIC:A9392-2"),
            ("BIOUG<CAN>:BP2010-044", "BIOUG<CAN>:BP2010-045"),
            ("KIOST-MB-CR255-4", "KIOST-MB-CR255-5"),
            ("CE438", "CE438_2"),
        ],
    )
    def test_different_suffix_not_merged(self, va, vb):
        kw = dict(Ref1Title=PAPER, collection_date="12-Jun-2019", lat_lon="66.55 N 33.10 E")
        df = pd.DataFrame([rec("A", "P c", "coi", va, "k1", **kw), rec("B", "P c", "16s", vb, "k2", **kw)])
        out, report = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 2
        assert report["decision"].str.contains("sub-identifier").any()

    @pytest.mark.parametrize(
        "va,vb",
        [
            ("COI_ZMMU_MSU_WS399", "WS399"),
            ("WS399_COI", "28S_ZMMU_WS399"),
            ("AM W.51461", "AM:W.51461"),
            ("personal collection:Christer Erseus:CE14346", "CE14346"),
            ("28S_ZMMU_MSU_WS16799_Pr21", "COI_ZMMU_MSU_WS16799"),  # Pr21: an extra sample code, not a lot number
        ],
    )
    def test_same_suffix_still_merged(self, va, vb):
        df = pd.DataFrame(
            [rec("A", "P c", "coi", va, "k1", Ref1Title=PAPER), rec("B", "P c", "16s", vb, "k2", Ref1Title=PAPER)]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1

    def test_candidate_from_isolate_when_voucher_differs(self):
        """GU901728: voucher SMNH<SWE>:108431, isolate CE438; PP994571: isolate CE438 only."""
        df = pd.DataFrame(
            [
                rec("GU901728", "H g", "12s", "SMNH<SWE>:108431", "k1", isolate="CE438", Ref1Title=PAPER),
                rec("PP994571", "H g", "28s", "", "k2", isolate="CE438", Ref1Title=PAPER),
            ]
        )
        out, report = reconcile_dataframe(df)
        assert len(report) == 1
        assert out["species_voucher_new"].nunique() == 1

    def test_merge_when_top_field_has_no_core(self):
        """voucher 's.n.' (no identifier) + isolate AB1234 meets isolate AB1234 (crashed on Nymphaeaceae)."""
        df = pd.DataFrame(
            [
                rec("A", "P c", "coi", "s.n.", "k1", isolate="AB1234", Ref1Title=PAPER),
                rec("B", "P c", "16s", "no voucher", "k2", isolate="AB1234", Ref1Title=PAPER),
            ]
        )
        out, _ = reconcile_dataframe(df)
        assert out["species_voucher_new"].nunique() == 1
        assert out["species_voucher_new"].iloc[0].endswith("AB1234")
