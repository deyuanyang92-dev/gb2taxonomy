"""Regression tests for the voucher.py / organize.py review fixes (L-08, L-09, L-11, L-12, L-22, L-23)."""

import importlib
import logging

import pandas as pd
import pytest

from g2t.organize import (
    OrganizeConfig,
    collect_locusids_per_gene,
    metadata_first_nonempty,
    organize,
    process_group,
)
from g2t.organize import main as organize_main
from g2t.voucher import build_species_vouchers
from g2t.voucher import main as voucher_main

# g2t.voucher is shadowed by a function of the same name in g2t/__init__; import the module explicitly.
voucher_mod = importlib.import_module("g2t.voucher")


def _read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


# ---------------------------------------------------------------------------
# L-08: mtgenome / 18-28s propagation honours config.gene_includes
# ---------------------------------------------------------------------------


class TestL08MtgenomeIncludes:
    def test_custom_mtgenome_includes_limits_propagation(self):
        group = pd.DataFrame({"gene_type": ["mtgenome"], "LocusID": ["MG1"]})
        cfg = OrganizeConfig(gene_includes={
            "18-28s": ["18s", "28s", "its1-its2"],
            "mtgenome": ["coi"],
        })
        res = collect_locusids_per_gene(group, cfg)
        assert res["mtgenome"] == "MG1"
        assert res["coi"] == "MG1"
        for g in ("16s", "12s", "cob", "cox2", "cox3"):
            assert res[g] == "", g

    def test_cli_mtgenome_includes_takes_effect(self, tmp_path):
        inp = tmp_path / "in.csv"
        inp.write_text(
            "species_voucher_new,organism,gene_type,LocusID\n"
            "Sp_A_1,Sp A,mtgenome,MG1\n"
        )
        out_default = tmp_path / "default.csv"
        out_coi = tmp_path / "coi_only.csv"
        organize_main(["-i", str(inp), "-o", str(out_default)])
        organize_main(["-i", str(inp), "-o", str(out_coi), "--mtgenome_includes", "coi"])
        d = _read(out_default).iloc[0]
        c = _read(out_coi).iloc[0]
        assert [d[g] for g in ("coi", "16s", "12s", "cob", "cox2", "cox3")] == ["MG1"] * 6
        assert c["coi"] == "MG1"
        assert [c[g] for g in ("16s", "12s", "cob", "cox2", "cox3")] == [""] * 5

    def test_default_propagation_unchanged(self):
        """Default config: mtgenome goes to every mito column of gene_order (get_mito_genes)."""
        group = pd.DataFrame({"gene_type": ["mtgenome"], "LocusID": ["MG1"]})
        res = collect_locusids_per_gene(group, OrganizeConfig())
        for g in ("coi", "16s", "12s", "cob", "cox2", "cox3"):
            assert res[g] == "MG1", g
        for g in ("18s", "28s", "its1-its2", "18-28s", "ef-1", "h3"):
            assert res[g] == "", g

    def test_default_follows_gene_order_mito_genes(self):
        """Pin existing default behaviour: a custom gene_order with nd1 still gets mtgenome copied to nd1."""
        group = pd.DataFrame({"gene_type": ["mtgenome"], "LocusID": ["MG1"]})
        cfg = OrganizeConfig(gene_order=["mtgenome", "coi", "nd1"])
        res = collect_locusids_per_gene(group, cfg)
        assert res["coi"] == "MG1"
        assert res["nd1"] == "MG1"

    def test_18_28s_uses_gene_includes(self):
        group = pd.DataFrame({"gene_type": ["18-28s"], "LocusID": ["R1"]})
        cfg = OrganizeConfig(gene_includes={
            "18-28s": ["18s"],
            "mtgenome": ["coi", "16s", "12s", "cob", "cox2", "cox3"],
        })
        res = collect_locusids_per_gene(group, cfg)
        assert res["18s"] == "R1"
        assert res["28s"] == ""
        assert res["its1-its2"] == ""
        default = collect_locusids_per_gene(group, OrganizeConfig())
        assert default["18s"] == default["28s"] == default["its1-its2"] == "R1"


# ---------------------------------------------------------------------------
# L-09: Conflict = any(); reason fields keep every record's value (with LocusID)
# ---------------------------------------------------------------------------


@pytest.fixture
def conflict_group():
    return pd.DataFrame({
        "species_voucher_new": ["Sp_A_1", "Sp_A_1"],
        "organism": ["Sp A", "Sp A"],
        "gene_type": ["coi", "16s,cox2"],
        "LocusID": ["A1.1", "A2.1"],
        "Conflict": ["False", "True"],
        "Conflict_reason": ["", "Multiple candidate gene types: 16s, cox2"],
        "Original_match": ["coi", "16s,cox2"],
        "Assignment_reason": ["coi: keyword", "16s: x; cox2: y"],
        "geo_loc_name": ["Russia: White Sea", ""],
        "country": ["", "Norway"],
    })


class TestL09ConflictAggregation:
    def test_conflict_is_any(self, conflict_group):
        res = metadata_first_nonempty(conflict_group, ["Conflict"])
        assert res["Conflict"] == "True"

    def test_conflict_true_first_still_true(self, conflict_group):
        res = metadata_first_nonempty(conflict_group.iloc[::-1], ["Conflict"])
        assert res["Conflict"] == "True"

    def test_conflict_all_false_stays_false(self, conflict_group):
        conflict_group["Conflict"] = ["False", "False"]
        assert metadata_first_nonempty(conflict_group, ["Conflict"])["Conflict"] == "False"

    def test_conflict_all_empty_stays_empty(self, conflict_group):
        conflict_group["Conflict"] = ["", ""]
        assert metadata_first_nonempty(conflict_group, ["Conflict"])["Conflict"] == ""

    def test_reason_fields_joined_with_locusid(self, conflict_group):
        res = metadata_first_nonempty(
            conflict_group, ["Conflict_reason", "Original_match", "Assignment_reason"])
        assert res["Conflict_reason"] == "A2.1: Multiple candidate gene types: 16s, cox2"
        assert res["Original_match"] == "A1.1: coi | A2.1: 16s,cox2"
        assert res["Assignment_reason"] == "A1.1: coi: keyword | A2.1: 16s: x; cox2: y"

    def test_identical_values_are_deduplicated(self, conflict_group):
        conflict_group["Original_match"] = ["coi", "coi"]
        res = metadata_first_nonempty(conflict_group, ["Original_match"])
        assert res["Original_match"].count("coi") == 1
        assert res["Original_match"].startswith("A1.1")

    def test_other_metadata_keep_first_nonempty(self, conflict_group):
        res = metadata_first_nonempty(conflict_group, ["geo_loc_name", "country"])
        assert res == {"geo_loc_name": "Russia: White Sea", "country": "Norway"}

    def test_process_group_and_organize_end_to_end(self, conflict_group, tmp_path):
        row = process_group(conflict_group, OrganizeConfig())
        assert row["Conflict"] == "True"
        assert "A2.1" in row["Conflict_reason"]
        inp = tmp_path / "in.csv"
        conflict_group.to_csv(inp, index=False)
        out = tmp_path / "out.csv"
        assert organize(str(inp), str(out)).success
        df = _read(out)
        assert df.loc[0, "Conflict"] == "True"
        assert df.loc[0, "coi"] == "A1.1"
        assert df.loc[0, "16s"] == "A2.1"

    def test_without_locusid_column_does_not_crash(self):
        g = pd.DataFrame({"Conflict_reason": ["x", "y"]})
        assert metadata_first_nonempty(g, ["Conflict_reason"])["Conflict_reason"] == "x | y"


# ---------------------------------------------------------------------------
# L-11: normalize_column_names only touches spaces, never the case
# ---------------------------------------------------------------------------


class TestL11NormalizeColumns:
    def test_case_preserved_spaces_replaced(self, tmp_path):
        inp = tmp_path / "in.csv"
        inp.write_text(
            "LocusID,ACCESSION,Organism,organism,gene type, Sample Id ,gene_type\n"
            "L1,A1,Org X,org x,coi,S1,coi\n"
        )
        out_dir = tmp_path / "out"
        res = build_species_vouchers(str(inp), str(out_dir), normalize_column_names=True, quiet=True)
        assert res.success
        cols = list(_read(out_dir / "updated_species_voucher.csv").columns)
        for c in ("LocusID", "ACCESSION", "Organism", "organism", "Sample_Id"):
            assert c in cols, c
        assert "locusid" not in cols
        assert "gene_type" in cols
        # "gene type" collides with "gene_type" -> deduplicated, never dropped
        assert any(c.startswith("gene_type__dup") for c in cols)

    def test_pipeline_with_normalize_columns_reaches_organize(self, tmp_path):
        import g2t

        out_root = tmp_path / "run"
        classify_dir = out_root / "labeled_genes"
        classify_dir.mkdir(parents=True)
        (out_root / "gb_metadata").mkdir(parents=True)
        (out_root / "gb_metadata" / "final.csv").write_text("LocusID\nA.1\n")
        pd.DataFrame({
            "gene_type": ["coi", "16s"],
            "Conflict": ["False", "False"],
            "LocusID": ["A.1", "B.1"],
            "ACCESSION": ["A", "B"],
            "Organism": ["Priapulus caudatus", "Priapulus caudatus"],
            "organism": ["Priapulus caudatus", "Priapulus caudatus"],
            "specimen voucher": ["V 1", "V 1"],
        }).to_csv(classify_dir / "assigned_genes_types_all.csv", index=False)
        # "specimen voucher" -> "specimen_voucher" must still feed the voucher key
        result = g2t.run(
            input_files=[], output_dir=str(out_root), quiet=True, normalize_columns=True,
            skip_extract=True, skip_classify=True, skip_reconcile=True,
        )
        assert result.success, result.log
        org = _read(out_root / "organized_genes" / "organized_species_voucher.csv")
        assert len(org) == 1
        assert org.loc[0, "species_voucher_new"] == "Priapulus_caudatus_V_1"
        assert org.loc[0, "coi"] == "A.1"
        assert org.loc[0, "16s"] == "B.1"


# ---------------------------------------------------------------------------
# L-12: CLI defaults == API defaults
# ---------------------------------------------------------------------------


@pytest.fixture
def defaults_csv(tmp_path):
    df = pd.DataFrame([
        dict(LocusID="A.1", ACCESSION="A", organism="Priapulus caudatus", gene_type="coi",
             culture_collection="CC-7", haplotype="H1"),
        dict(LocusID="B.1", ACCESSION="B", organism="Priapulus caudatus", gene_type="16s",
             culture_collection="CC-7", haplotype="H1"),
        dict(LocusID="C.1", ACCESSION="C", organism="Priapulus caudatus", gene_type="coi",
             haplotype="H9"),
        dict(LocusID="D.1", ACCESSION="D", organism="Priapulus caudatus", gene_type="coi"),
    ])
    p = tmp_path / "in.csv"
    df.to_csv(p, index=False)
    return p


class TestL12CliDefaults:
    def test_default_voucher_columns_constant(self):
        cols = voucher_mod.DEFAULT_VOUCHER_COLUMNS
        assert list(cols) == ["specimen_voucher", "isolate", "culture_collection", "clone", "strain"]

    def test_cli_matches_api(self, defaults_csv, tmp_path):
        api_dir, cli_dir = tmp_path / "api", tmp_path / "cli"
        build_species_vouchers(str(defaults_csv), str(api_dir), quiet=True)
        voucher_main(["-i", str(defaults_csv), "-o", str(cli_dir), "--quiet"])
        api = _read(api_dir / "updated_species_voucher.csv")
        cli = _read(cli_dir / "updated_species_voucher.csv")
        assert list(api["species_voucher_new"]) == list(cli["species_voucher_new"])
        assert list(cli["species_voucher_new"]) == [
            "Priapulus_caudatus_CC-7", "Priapulus_caudatus_CC-7",
            "Priapulus_caudatus_C", "Priapulus_caudatus_D"]

    def test_cli_fill_haplotype_still_opt_in(self, defaults_csv, tmp_path):
        out = tmp_path / "hap"
        voucher_main(["-i", str(defaults_csv), "-o", str(out), "--quiet", "--fill_haplotype", "true"])
        new = list(_read(out / "updated_species_voucher.csv")["species_voucher_new"])
        assert new[2] == "Priapulus_caudatus_H9"

    def test_organize_cli_metadata_columns_default_is_config_default(self, tmp_path):
        inp = tmp_path / "in.csv"
        inp.write_text(
            "species_voucher_new,organism,gene_type,LocusID,country,collection_date,isolate,"
            "specimen_voucher,host,habitat,isolation_source,culture_collection,clone,strain,note\n"
            "Sp_A_1,Sp A,coi,A1.1,Norway,2001,I1,V1,H,Hab,Src,CC,Cl,St,Nt\n"
        )
        out = tmp_path / "out.csv"
        organize_main(["-i", str(inp), "-o", str(out)])
        df = _read(out)
        expected = OrganizeConfig().meta_columns
        for c in expected:
            assert c in df.columns, c
        assert df.loc[0, "country"] == "Norway"
        assert df.loc[0, "collection_date"] == "2001"
        assert df.loc[0, "strain"] == "St"


# ---------------------------------------------------------------------------
# L-22: whitespace / punctuation-only voucher values count as empty
# ---------------------------------------------------------------------------


class TestL22BlankVoucher:
    @pytest.fixture
    def edge_csv(self, tmp_path):
        rows = [
            dict(LocusID="A.1", ACCESSION="AA", organism="Priapulus caudatus", gene_type="coi",
                 specimen_voucher=" ", isolate="I1"),
            dict(LocusID="B.1", ACCESSION="BB", organism="Priapulus caudatus", gene_type="coi",
                 specimen_voucher="..", isolate="I2"),
            dict(LocusID="C.1", ACCESSION="CC", organism="Priapulus caudatus", gene_type="16s",
                 specimen_voucher="()", isolate="I3"),
            dict(LocusID="D.1", ACCESSION="DD", organism="Priapulus caudatus", gene_type="coi",
                 specimen_voucher="-", isolate=" . "),
            dict(LocusID="E.1", ACCESSION="EE", organism="Priapulus caudatus", gene_type="coi",
                 specimen_voucher="ab 1", isolate="I5"),
        ]
        p = tmp_path / "edge.csv"
        pd.DataFrame(rows).to_csv(p, index=False, na_rep="")
        return p

    def test_blank_values_fall_back_to_next_column_then_accession(self, edge_csv, tmp_path):
        out = tmp_path / "out"
        build_species_vouchers(str(edge_csv), str(out), quiet=True)
        df = _read(out / "updated_species_voucher.csv")
        assert list(df["species_voucher_new"]) == [
            "Priapulus_caudatus_I1",
            "Priapulus_caudatus_I2",
            "Priapulus_caudatus_I3",
            "Priapulus_caudatus_DD",
            "Priapulus_caudatus_ab_1",
        ]

    def test_pse_keeps_raw_value_for_valid_entries(self, edge_csv, tmp_path):
        out = tmp_path / "out"
        build_species_vouchers(str(edge_csv), str(out), quiet=True)
        df = _read(out / "updated_species_voucher.csv")
        assert df.loc[4, "species_voucher_pse"] == "ab 1"
        assert df.loc[0, "species_voucher_pse"] == "I1"

    def test_blank_haplotype_does_not_block_accession_fallback(self, tmp_path):
        p = tmp_path / "hap.csv"
        pd.DataFrame([dict(LocusID="A.1", ACCESSION="AA", organism="Sp x", gene_type="coi",
                           haplotype="..")]).to_csv(p, index=False)
        out = tmp_path / "out"
        build_species_vouchers(str(p), str(out), if_write_haplotype=True, quiet=True)
        assert _read(out / "updated_species_voucher.csv").loc[0, "species_voucher_new"] == "Sp_x_AA"


# ---------------------------------------------------------------------------
# L-23: dropped gene types are reported; empty group keys are not merged
# ---------------------------------------------------------------------------


class TestL23UnknownGenesAndKeys:
    def test_dropped_gene_types_are_warned_with_counts(self, tmp_path, caplog):
        inp = tmp_path / "in.csv"
        inp.write_text(
            "species_voucher_new,organism,gene_type,LocusID\n"
            "Sp_A_1,Sp A,nd1,N1.1\n"
            "Sp_B_1,Sp B,\"nd1,coi\",N2.1\n"
            "Sp_C_1,Sp C,nd6#,N3.1\n"
            "Sp_D_1,Sp D,coi,C1.1\n"
        )
        out = tmp_path / "out.csv"
        with caplog.at_level(logging.WARNING, logger="g2t.organize"):
            assert organize(str(inp), str(out)).success
        msgs = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
        joined = " ".join(msgs)
        assert "nd1=2" in joined
        assert "nd6=1" in joined
        assert "coi=" not in joined

    def test_no_warning_when_all_genes_known(self, tmp_path, caplog):
        inp = tmp_path / "in.csv"
        inp.write_text(
            "species_voucher_new,organism,gene_type,LocusID\n"
            "Sp_A_1,Sp A,coi,A1.1\n"
            "Sp_A_1,Sp A,\"16s,18s\",A2.1\n"
        )
        with caplog.at_level(logging.WARNING, logger="g2t.organize"):
            organize(str(inp), str(tmp_path / "out.csv"))
        assert not [r for r in caplog.records if r.levelno >= logging.WARNING]

    def test_custom_gene_in_gene_order_is_kept_and_not_warned(self, tmp_path, caplog):
        inp = tmp_path / "in.csv"
        inp.write_text(
            "species_voucher_new,organism,gene_type,LocusID\n"
            "Sp_A_1,Sp A,nd1,N1.1\n"
        )
        cfg = OrganizeConfig(gene_order=["coi", "nd1"])
        with caplog.at_level(logging.WARNING, logger="g2t.organize"):
            organize(str(inp), str(tmp_path / "out.csv"), cfg)
        assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
        assert _read(tmp_path / "out.csv").loc[0, "nd1"] == "N1.1"

    def test_empty_group_keys_become_separate_unknown_rows(self, tmp_path):
        inp = tmp_path / "in.csv"
        inp.write_text(
            "species_voucher_new,organism,gene_type,LocusID\n"
            ",Sp X,coi,X1.1\n"
            ",Sp Y,16s,Y1.1\n"
            "Sp_D_1,Sp D,coi,D1.1\n"
        )
        out = tmp_path / "out.csv"
        assert organize(str(inp), str(out)).success
        df = _read(out).set_index("species_voucher_new")
        assert len(df) == 3
        assert "UNKNOWN" not in df.index
        assert df.loc["UNKNOWN_X1.1", "organism"] == "Sp X"
        assert df.loc["UNKNOWN_X1.1", "coi"] == "X1.1"
        assert df.loc["UNKNOWN_X1.1", "16s"] == ""
        assert df.loc["UNKNOWN_Y1.1", "organism"] == "Sp Y"
        assert df.loc["UNKNOWN_Y1.1", "16s"] == "Y1.1"
