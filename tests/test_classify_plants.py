"""Plant (plastid) gene types; DEFINITION lines taken from Nymphaeaceae, Ginkgoaceae and Ophioglossaceae records."""

import pytest

from g2t.classify import MatchConfig, process_row


@pytest.fixture
def config():
    return MatchConfig()


def typed(config, definition, length, organelle="plastid:chloroplast"):
    row = {"Definition": definition, "organelle": organelle, "Length": str(length), "LocusID": "X1"}
    return process_row(row, config)["gene_type"]


@pytest.mark.parametrize(
    "definition,length,expected",
    [
        (
            "Nymphaea alba ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit (rbcL) gene, partial cds; "
            "chloroplast",
            1350,
            "rbcl",
        ),
        ("Nymphaea alba maturase K (matK) gene, partial cds; chloroplast", 800, "matk"),
        ("Nuphar lutea maturase K-like (matK) gene, partial sequence; chloroplast", 750, "matk"),
        ("Nuphar lutea trnL-trnF intergenic spacer, partial sequence; chloroplast", 400, "trnl-trnf"),
        (
            "Nuphar lutea trnT-trnL intergenic spacer, partial sequence; tRNA-Leu (trnL) gene, complete sequence; and "
            "trnL-trnF intergenic spacer, partial sequence; chloroplast",
            1500,
            "trnl-trnf",
        ),
        (
            "Nuphar lutea PsbA (psbA) gene, partial cds; and psbA-trnH intergenic spacer, partial sequence; "
            "chloroplast",
            450,
            "psba-trnh",
        ),
        (
            "Nuphar lutea tRNA-His (trnH) gene, partial sequence; trnH-psbA intergenic spacer, complete sequence; and "
            "PsbA (psbA) gene, partial cds; chloroplast",
            500,
            "psba-trnh",
        ),
        ("Nuphar lutea ribosomal protein L16 (rpl16) gene, intron; chloroplast", 1000, "rpl16"),
        ("Nuphar lutea rpL16 gene, intron; chloroplast", 1000, "rpl16"),
        ("Ginkgo biloba chloroplast, complete genome", 156945, "plastome"),
        ("Ophioglossum vulgatum plastid, complete genome", 138270, "plastome"),
    ],
)
def test_plastid_markers(config, definition, length, expected):
    assert typed(config, definition, length) == expected


def test_spacer_name_is_not_the_gene(config):
    """'rbcL-atpB intergenic spacer' names a spacer (type rbcl-atpb), not the rbcL gene."""
    d = (
        "Nuphar lutea rbcL-atpB intergenic spacer, partial sequence; and ATP synthase beta subunit (atpB) gene, "
        "partial cds; chloroplast"
    )
    assert typed(config, d, 900) == "rbcl-atpb"


def test_plastid_word_in_definition_without_organelle(config):
    d = "Nymphaea alba ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit (rbcL) gene, partial cds; plastid"
    assert typed(config, d, 1300, organelle="") == "rbcl"


def test_nuclear_gene_for_chloroplast_product_is_nuclear(config):
    d = (
        "Nymphaea alba glyceraldehyde-3-phosphate dehydrogenase (GapA) gene, partial cds; "
        "nuclear gene for chloroplast product"
    )
    assert typed(config, d, 900, organelle="") == ""


def test_plant_nuclear_26s_is_28s(config):
    assert typed(config, "Nymphaea alba 26S ribosomal RNA gene, partial sequence", 1200, organelle="") == "28s"


def test_plant_its_unchanged(config):
    d = (
        "Nymphaea alba internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene, complete sequence; "
        "and internal transcribed spacer 2, partial sequence"
    )
    assert typed(config, d, 700, organelle="") == "its1-its2"


def test_plastome_not_removed_by_length_filter(tmp_path):
    from g2t.classify import classify

    inp = tmp_path / "final.csv"
    inp.write_text(
        "LocusID,ACCESSION,Definition,Length,organelle,MoleculeType,Topology\n"
        "NC1.1,NC1,Ginkgo biloba chloroplast complete genome,156945 bp,plastid:chloroplast,DNA,circular\n"
    )
    res = classify(str(inp), str(tmp_path / "out"))
    assert res.success
    import pandas as pd

    st = pd.read_csv(tmp_path / "out" / "record_status.csv")
    assert st.loc[0, "status"] == "assigned" and st.loc[0, "gene_type"] == "plastome"


def test_matrix_columns_and_plastome_propagation():
    import pandas as pd

    from g2t.organize import OrganizeConfig, collect_locusids_per_gene

    cfg = OrganizeConfig()
    for g in ("plastome", "rbcl", "matk", "trnl-trnf", "psba-trnh"):
        assert g in cfg.gene_order
    grp = pd.DataFrame({"gene_type": ["plastome", "its1-its2"], "LocusID": ["NC1.1", "ITS1.1"]})
    ids = collect_locusids_per_gene(grp, cfg)
    assert ids["rbcl"] == "NC1.1" and ids["matk"] == "NC1.1"
    assert ids["coi"] == ""


@pytest.mark.parametrize(
    "definition,expected",
    [
        (
            "Nuphar lutea ribosomal protein L32 (rpl32) gene, partial cds; and rpl32-trnL intergenic spacer, partial "
            "sequence; chloroplast",
            "rpl32-trnl",
        ),
        ("Ginkgo biloba trnS-trnG intergenic spacer, partial sequence; chloroplast", "trns-trng"),
        (
            "Ophioglossum sp. AtpF (atpF) gene, partial cds; atpF-atpH intergenic spacer, complete sequence; and AtpH "
            "(atpH) gene, partial cds; chloroplast",
            "atpf-atph",
        ),
        ("Nuphar lutea petB-petD intergenic spacer region, partial sequence; plastid", "petb-petd"),
        (
            "Nuphar lutea PsbE (psbE) gene, partial cds; psbE-petL intergenic spacer, complete sequence; "
            "and PetL (petL) gene, partial cds; chloroplast",
            "psbe-petl",
        ),
        (
            "Ophioglossum sp. rbcL-atpB intergenic spacer, partial sequence; "
            "and ATP synthase beta subunit (atpB) gene, partial cds; chloroplast",
            "rbcl-atpb",
        ),
        ("Nuphar lutea trnK gene, introns; chloroplast", "trnk"),
        (
            "Nuphar lutea tRNA-Lys (trnK) gene, partial sequence; "
            "and maturase K (matK) gene, complete cds; chloroplast",
            "matk",
        ),
    ],
)
def test_plastid_spacers(config, definition, expected):
    assert typed(config, definition, 900) == expected


# --- review round 1 (reviewer 2, M5) ---


def test_plant_mito_fragment_is_not_a_mitogenome(config):
    """Ginkgo MW354248: 5.2 kb with rrn26 + trnY; 'complete sequence' alone does not make a mitogenome."""
    d = "Ginkgo biloba 26S ribosomal RNA (rrn26) and tRNA-Tyr (trnY) genes, complete sequence; mitochondrial"
    gt = typed(config, d, 5224, organelle="mitochondrion")
    assert "mtgenome" not in gt.split(",")


def test_mito_rrn26_not_put_in_nuclear_28s():
    """A mitochondrial 26S rRNA must not land in the nuclear 28s column via round 2."""
    from g2t.classify import recheck_match_row

    row = {
        "Definition": "Ginkgo biloba 26S ribosomal RNA (rrn26) gene, partial sequence; mitochondrial",
        "organelle": "mitochondrion",
        "Length": "1500",
        "MoleculeType": "DNA",
        "Topology": "linear",
        "Taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Ginkgoopsida",
    }
    assert "28s" not in recheck_match_row(row, MatchConfig())["gene_type"].split(",")


def test_animal_18s_titled_mitochondrial_keeps_gene_and_is_flagged():
    """Phyllodocida KT270941 '18S ribosomal RNA gene, partial sequence; mitochondrial': a submitter's slip."""
    import tempfile
    from pathlib import Path

    import pandas as pd

    from g2t.classify import classify

    d = Path(tempfile.mkdtemp())
    (d / "final.csv").write_text(
        "LocusID,ACCESSION,Definition,Length,organelle,MoleculeType,Topology,Taxonomy\n"
        'KT270941.1,KT270941,"Namalycastis abiuma isolate TRL423 18S ribosomal RNA gene, partial sequence; '
        'mitochondrial",1700 bp,mitochondrion,DNA,linear,Eukaryota; Metazoa; Annelida\n'
    )
    classify(str(d / "final.csv"), str(d / "out"))
    st = pd.read_csv(d / "out" / "record_status.csv")
    assert st.loc[0, "gene_type"] == "18s"


@pytest.mark.parametrize(
    "definition,organelle",
    [
        ("Ginkgo biloba chloroplast DNA, complete genome", "plastid:chloroplast"),
        ("Nuphar lutea genome assembly, organelle: plastid:chloroplast", "plastid:chloroplast"),
    ],
)
def test_plastome_titles(config, definition, organelle):
    assert typed(config, definition, 156945, organelle=organelle) == "plastome"


def test_circular_plastid_is_plastome_in_round_2():
    from g2t.classify import recheck_match_row

    row = {
        "Definition": "Ginkgo biloba chloroplast DNA, complete genome",
        "organelle": "plastid:chloroplast",
        "Length": "156945",
        "MoleculeType": "DNA",
        "Topology": "circular",
    }
    assert recheck_match_row(row, MatchConfig())["gene_type"] == "plastome"


@pytest.mark.parametrize(
    "definition,expected",
    [
        (
            "Botrychium lunaria ribosomal protein S4 (rps4) gene, partial cds; and rps4-trnS intergenic spacer, "
            "partial sequence; chloroplast",
            "rps4",
        ),
        ("Ophioglossum sp. ATP synthase beta subunit (atpB) gene, partial cds; chloroplast", "atpb"),
    ],
)
def test_fern_markers(config, definition, expected):
    assert typed(config, definition, 900) == expected


def test_animal_mitogenome_titles_still_typed(config):
    for d in (
        "Priapulus caudatus mitochondrion, complete genome",
        "Hesione sp. mitochondrial DNA, complete genome",
        "UNVERIFIED: Hesione genetta mitochondrion sequence",
    ):
        assert typed(config, d, 15500, organelle="mitochondrion") == "mtgenome", d
