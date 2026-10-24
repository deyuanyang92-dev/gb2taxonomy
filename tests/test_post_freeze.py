"""Problems found after the rules were frozen (held-out tables and blinded audits); examples are real records."""

import pandas as pd
import pytest

from g2t.classify import MatchConfig, process_row
from g2t.reconcile import reconcile_dataframe
from g2t.voucher import _is_blank_value


def typed(definition, length, organelle=""):
    row = {"Definition": definition, "organelle": organelle, "Length": str(length), "LocusID": "X1"}
    return set(process_row(row, MatchConfig())["gene_type"].split(",")) - {""}


@pytest.mark.parametrize(
    "definition,length,organelle,expected",
    [
        # 30: a spacer record that also names a gene (Angraecum KF672272)
        (
            "Angraecum cf. humblotianum trnL gene, partial sequence; and maturase K (matK) gene, partial cds; "
            "chloroplast",
            1763,
            "plastid:chloroplast",
            {"trnl-trnf", "matk"},
        ),
        # 34: unannotated plastomes (Leptochilus MW876331)
        ("UNVERIFIED: Leptochilus decurrens chloroplast sequence", 156868, "plastid:chloroplast", {"plastome"}),
        ("UNVERIFIED: Leptochilus longissimus plastid sequence", 156375, "plastid", {"plastome"}),
        # 37: 'gene for LSU/SSU rRNA' (Opheliidae ICHUM records)
        ("Euzonus ezoensis ICHUM:9159 gene for LSU rRNA, partial sequence", 800, "", {"28s"}),
        ("Euzonus ezoensis ICHUM:9168 gene for SSU rRNA, partial sequence", 1700, "", {"18s"}),
        (
            "Trypanedenta gigantea isolate Tgi_S20567 18S small ribosomal RNA subunit gene, partial sequence",
            1772,
            "",
            {"18s"},
        ),
        # 38: plant 26S flank of an ITS record (Campylocentrum micranthum)
        (
            "Campylocentrum micranthum isolate 2 18S ribosomal RNA gene, partial sequence; internal transcribed spacer "
            "1, 5.8S ribosomal RNA gene and internal transcribed spacer 2, complete sequence; and 26S ribosomal RNA "
            "gene, partial sequence",
            677,
            "",
            {"18s", "its1-its2", "28s"},
        ),
    ],
)
def test_gene_types(definition, length, organelle, expected):
    assert typed(definition, length, organelle) == expected


def test_short_chloroplast_sequence_is_not_a_plastome():
    assert "plastome" not in typed("UNVERIFIED_ASMBLY: Ophioglossum sp. chloroplast sequence", 286, "plastid")


@pytest.mark.parametrize("v", ["BRLU", "KUN", "MNHN"])  # 35: an institution code alone identifies no specimen
def test_voucher_without_number_is_blank(v):
    assert _is_blank_value(v, specimen_voucher=True)
    assert not _is_blank_value(v)  # as an isolate or haplotype label it is kept


def rec(acc, org, gene, voucher, key, **kw):
    row = dict(
        ACCESSION=acc,
        organism=org,
        gene_type=gene,
        specimen_voucher=voucher,
        isolate="",
        species_voucher_new=key,
        Ref1Title="Same paper",
        Ref2Title="",
        Ref1Authors="",
        collection_date="12-Nov-2011",
        lat_lon="",
        geo_loc_name="",
        country="",
        collected_by="",
    )
    row.update(kw)
    return row


def test_different_specimen_numbers_around_a_shared_code():
    """36: Polyeunoa laevis NSJC66_1005_S001 / NSJC66_962_S001 are two specimens."""
    df = pd.DataFrame(
        [
            rec("A", "Polyeunoa laevis", "coi", "", "k1", isolate="NSJC66_1005_S001"),
            rec("B", "Polyeunoa laevis", "16s", "", "k2", isolate="NSJC66_962_S001"),
        ]
    )
    out, _ = reconcile_dataframe(df)
    assert out["species_voucher_new"].nunique() == 2


def test_individuals_of_a_split_lot_are_not_merged_back():
    """40: MCZ 70173 isolates 5 and 7 (different genes) were split as a lot, then merged on the identical voucher."""
    rows = [
        rec("C5", "B s", "coi", "MCZ 70173", "k", isolate="5"),
        rec("C6", "B s", "coi", "MCZ 70173", "k", isolate="6"),
        rec("S7", "B s", "16s", "MCZ 70173", "k", isolate="7"),
    ]
    out, _ = reconcile_dataframe(pd.DataFrame(rows))
    assert out["species_voucher_new"].nunique() == 3


def test_ignore_evidence_for_sensitivity_analysis():
    from g2t.reconcile import ReconcileConfig

    df = pd.DataFrame(
        [
            rec("A", "P c", "coi", "COI_ZMMU_WS399", "k1", collection_date=""),
            rec("B", "P c", "28s", "28S_ZMMU_WS399", "k2", collection_date=""),
        ]
    )
    assert reconcile_dataframe(df)[0]["species_voucher_new"].nunique() == 1  # shared paper
    cfg = ReconcileConfig(ignore_evidence=("same publication",))
    assert reconcile_dataframe(df, cfg)[0]["species_voucher_new"].nunique() == 2


def test_drop_empty_gene_columns(tmp_path):
    from g2t.organize import OrganizeConfig, organize

    inp = tmp_path / "in.csv"
    inp.write_text("species_voucher_new,organism,gene_type,LocusID\nk1,P c,coi,A.1\nk1,P c,16s,B.1\n")
    out = tmp_path / "m.csv"
    organize(str(inp), str(out))
    assert "rbcl" in pd.read_csv(out).columns  # default: every gene column, a fixed layout
    organize(str(inp), str(out), OrganizeConfig(drop_empty_gene_columns=True))
    cols = pd.read_csv(out).columns
    assert "coi" in cols and "16s" in cols and "rbcl" not in cols and "12s" not in cols


def test_pipeline_cli_has_drop_empty_option():
    from g2t._cli import main_pipeline

    with pytest.raises(SystemExit):
        main_pipeline(["--help"])
