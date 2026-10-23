"""Linking cases raised by the reviewers (round 1); examples are the accessions they cited."""

import pandas as pd
import pytest

from g2t.reconcile import reconcile_dataframe
from g2t.voucher import _is_blank_value

PAPER = "Same paper"


def rec(acc, org, gene, voucher, key, **kw):
    row = dict(
        ACCESSION=acc,
        organism=org,
        gene_type=gene,
        specimen_voucher=voucher,
        isolate="",
        species_voucher_new=key,
        Ref1Title=PAPER,
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


@pytest.mark.parametrize("v", ["FCO:n./d.", "n/d", "s.n.", "no voucher", "NA", "unvouchered", "unknown", "-"])
def test_placeholder_vouchers_are_blank(v):
    assert _is_blank_value(v)


@pytest.mark.parametrize("v", ["Juslen s.n. 2003", "SMNH 108431", "CE438", "FCO:1234"])
def test_real_vouchers_are_not_blank(v):
    assert not _is_blank_value(v)


@pytest.mark.parametrize(
    "va,vb",
    [
        ("USNM 123456", "MCZ 123456"),  # different institutions
        ("Christenhusz 5678 (K)", "Hassemer 5678 (FLOR)"),  # different collectors
        ("Juslen s.n. 2003", "Muola s.n. 2003"),  # different collectors; 2003 is a year (DQ026594 / DQ026627)
        ("Polychaeta1D9_St95", "Polychaeta1G8_St95"),  # two plate wells of one station (KF713428 / KF713468)
    ],
)
def test_different_institution_collector_or_sample_not_merged(va, vb):
    kw = dict(collection_date="12-Jun-2019")
    df = pd.DataFrame([rec("A", "P c", "coi", va, "k1", **kw), rec("B", "P c", "16s", vb, "k2", **kw)])
    out, report = reconcile_dataframe(df)
    assert out["species_voucher_new"].nunique() == 2, report["decision"].tolist()


def test_herbarium_collector_abbreviation_still_joins():
    """'C.B. Hellquist, J.H. Wiersema 16149 (NASC)' and 'H & W 16149; GH' (reviewer 2: correct join)."""
    df = pd.DataFrame(
        [
            rec("A", "Nuphar sp.", "matk", "C.B. Hellquist, J.H. Wiersema 16149 (NASC)", "k1"),
            rec("B", "Nuphar sp.", "rbcl", "H & W 16149; GH", "k2"),
        ]
    )
    out, _ = reconcile_dataframe(df)
    assert out["species_voucher_new"].nunique() == 1


def test_lot_with_several_isolates_of_one_gene_is_split():
    """MCZ 70173: 7 COI from different individuals of one lot (MH115399-MH115405)."""
    rows = [rec(f"MH11539{i}", "B s", "coi", "MCZ 70173", "Bs_MCZ_70173", isolate=f"ind{i}") for i in range(3)]
    rows.append(rec("X16S", "B s", "16s", "MCZ 70173", "Bs_MCZ_70173", isolate="ind0"))
    out, _ = reconcile_dataframe(pd.DataFrame(rows))
    key = out.set_index("ACCESSION")["species_voucher_new"]
    assert key.nunique() == 3
    assert key["X16S"] == key["MH115390"]  # 16S goes with the COI of the same isolate
    assert out["match_evidence"].str.contains("several individuals").all()


def test_one_voucher_one_isolate_per_gene_not_split():
    rows = [
        rec("A", "B s", "coi", "MCZ 70173", "k", isolate="ind1"),
        rec("B", "B s", "16s", "MCZ 70173", "k", isolate="ind2"),  # different isolate code, different gene
    ]
    out, _ = reconcile_dataframe(pd.DataFrame(rows))
    assert out["species_voucher_new"].nunique() == 1


def test_matrix_name_is_the_most_frequent():
    from g2t.organize import OrganizeConfig, process_group

    g = pd.DataFrame(
        {
            "species_voucher_new": ["k"] * 3,
            "organism": ["Botrychium pedunculosum", "Botrychium matricariifolium", "Botrychium matricariifolium"],
            "gene_type": ["matk", "rbcl", "trnl-trnf"],
            "LocusID": ["A.1", "B.1", "C.1"],
        }
    )
    out = process_group(g, OrganizeConfig())
    assert out["organism"] == "Botrychium matricariifolium"
    assert "Botrychium pedunculosum" in out["organism_names"]
