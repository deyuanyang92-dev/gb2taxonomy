"""g2t-mitostats: sequencing category, assembler names, publication year/journal (strings as found in GenBank)."""

import json
import re
from pathlib import Path

import pandas as pd
import pytest

from g2t import mitostats as ms

# ---------------------------------------------------------------- sequencing technology


@pytest.mark.parametrize(
    "raw, cat",
    [
        ("Illumina", "NGS"),
        ("HiSeqX Ten", "NGS"),
        ("DNBSEQTM platform (BGI-Shenzhen, China)", "NGS"),
        ("BGISEQ-500", "NGS"),
        ("MGI-2000", "NGS"),
        ("454", "NGS"),
        ("Ion Torrent", "NGS"),
        ("Sanger dideoxy sequencing", "Sanger"),
        ("ABI 3730xl", "Sanger"),
        ("PacBio", "Long-read"),
        ("Oxford Nanopore Technology", "Long-read"),
        ("Illumina; Nanopore", "Long-read"),  # hybrid assemblies count as long-read (thesis Fig. 2-12)
        ("Illumina; Sanger dideoxy sequencing", "NGS"),  # Sanger used to close gaps
        ("Sanger dideoxy sequencing; 454", "NGS"),
        ("", "Unknown"),
        ("unknown platform xyz", "Unknown"),
    ],
)
def test_sequencing_category(raw, cat):
    assert ms.sequencing_category(raw) == cat


def test_hybrid_is_flagged():
    assert ms.sequencing_platforms("Illumina; Oxford Nanopore") == ["Illumina", "Nanopore"]


# ---------------------------------------------------------------- assemblers


@pytest.mark.parametrize(
    "raw, names",
    [
        ("NOVOplasty v. 4.2.1", ["NOVOPlasty"]),
        ("mitoz v. v2.4; getorangelle v. v1.6.0; novoplasty v. v3.4.1", ["MitoZ", "GetOrganelle", "NOVOPlasty"]),
        ("GetOrgnelle v. 1.7.5.2", ["GetOrganelle"]),
        ("Trimmomatic v. 0.39; MitoFinder v. 1.4; MITOS v. 1", ["MitoFinder"]),
        ("SPADes v. 3.13.0; NOVOPlasty v. 4.3.1", ["SPAdes", "NOVOPlasty"]),
        ("IDBA_UD v. 03/15/2016; CLC Genomics Workbench v. 5.1", ["IDBA-UD", "CLC"]),
        ("GS newbler v. 2.9", ["Newbler"]),
        ("GS Data Analysis v. 2015", ["Newbler"]),
        ("Soap denovo v. 2", ["SOAPdenovo"]),
        ("novoplasty and mitofinder v. v4.3.1", ["NOVOPlasty", "MitoFinder"]),
        ("seqman v. v. Pro", ["SeqMan"]),
        ("DNASTAR v. 5.0", ["SeqMan"]),
        ("ABySS v. v2.0.2; GapCloser v. v1.12", ["ABySS"]),
        ("Getorganelle v. GetOrganelle v1.7.7.1", ["GetOrganelle"]),
        ("FastQC v. v4.2", []),
        ("", []),
    ],
)
def test_assemblers(raw, names):
    assert ms.assemblers(raw)[0] == names


def test_unrecognised_tool_is_kept():
    names, other = ms.assemblers("FooAssembler v. 1.2; NOVOPlasty v. 4.3.1")
    assert names == ["NOVOPlasty"] and other == ["FooAssembler v. 1.2"]
    assert ms.assemblers("2791862 v. v1.7.6.1")[1] == []  # no name, only numbers: not a program


@pytest.mark.parametrize(
    "names, other, raw, group",
    [
        (["NOVOPlasty"], [], "NOVOPlasty v. 4", "NOVOPlasty"),
        (["MitoZ", "NOVOPlasty"], [], "x", "Multi-assembly"),
        ([], ["weird"], "weird", "Other"),
        ([], [], "", "Not stated"),
        ([], [], "FastQC v. 1", "Not stated"),  # only non-assembly tools
    ],
)
def test_assembler_group(names, other, raw, group):
    assert ms.assembler_group(names, other) == group


# ---------------------------------------------------------------- references


@pytest.mark.parametrize(
    "journal, name, year",
    [
        ("Zool Stud 65, e37 (2026)", "Zool Stud", 2026),
        ("Mitochondrial DNA B Resour 5 (3), 2370-2371 (2020)", "Mitochondrial DNA B Resour", 2020),
        ("Mol. Phylogenet. Evol. 103, 1-10 (2016)", "Mol Phylogenet Evol", 2016),
        ("Genome Biol Evol (2020) In press", "Genome Biol Evol", 2020),
        ("Unpublished", "", None),
        ("Submitted (19-OCT-2023) Laboratory of Genetics", "", None),
    ],
)
def test_journal_name_and_year(journal, name, year):
    assert ms.journal_name_year(journal) == (name, year)


def test_publication_takes_first_published_ref():
    refs = [
        ("Unpublished", "T1"),
        ("Zool Stud 65, e37 (2026)", "T2"),
        ("Submitted (09-JUN-2023) Somewhere", ""),
    ]
    p = ms.publication(refs, create_date="2023-06-10")
    assert p == {"journal": "Zool Stud", "pub_status": "published", "pub_year": 2026, "submit_year": 2023,
                 "year": 2026, "year_source": "publication"}


def test_unpublished_falls_back_to_submission_year():
    p = ms.publication([("Unpublished", ""), ("Submitted (09-JUN-2023) X", "")], create_date="")
    assert p["pub_status"] == "unpublished" and p["year"] == 2023 and p["year_source"] == "submission"


# ---------------------------------------------------------------- end to end on header-only records

REC = """LOCUS       MN000001                   1 bp    DNA     linear   INV 01-JAN-2021
DEFINITION  Alitta virens mitochondrion, complete genome.
ACCESSION   MN000001
VERSION     MN000001.1
SOURCE      mitochondrion Alitta virens
  ORGANISM  Alitta virens
            Eukaryota; Metazoa; Spiralia; Lophotrochozoa; Annelida; Polychaeta;
            Errantia; Phyllodocida; Nereididae; Alitta.
REFERENCE   1  (bases 1 to 1)
  AUTHORS   Smith,A. and Lee,B.
  TITLE     The mitogenome of Alitta virens
  JOURNAL   Mitochondrial DNA B Resour 5 (3), 2370-2371 (2020)
REFERENCE   2  (bases 1 to 1)
  AUTHORS   Smith,A.
  TITLE     Direct Submission
  JOURNAL   Submitted (01-JAN-2020) Somewhere
COMMENT     ##Assembly-Data-START##
            Assembly Method       :: mitoz v. v2.4; Novoplasty v. v4.3.1
            Sequencing Technology :: Illumina
            ##Assembly-Data-END##
FEATURES             Location/Qualifiers
     source          1
                     /organism="Alitta virens"
                     /organelle="mitochondrion"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:6346"
ORIGIN
        1 a
//
"""


def test_table_from_header_records(tmp_path):
    d = tmp_path / "mitogenome-header"
    d.mkdir()
    (d / "batch_0001.gb").write_text(REC + REC.replace("MN000001", "MN000002").replace(
        "##Assembly-Data-START##\n            Assembly Method       :: mitoz v. v2.4; Novoplasty v. v4.3.1\n"
        "            Sequencing Technology :: Illumina\n            ##Assembly-Data-END##",
        "Sequenced by Sanger."))
    (d / "summary.tsv").write_text(
        "accession\tlength\tcreate_date\tupdate_date\ttaxid\ttitle\n"
        "MN000001.1\t15720\t2020/01/05\t2020/02/01\t6346\tx\n"
        "MN000002.1\t15710\t2020/01/05\t2020/02/01\t6346\tx\n".replace("/", "-")
    )
    lineage = {"6346": {"phylum": "Annelida", "class": "Polychaeta", "order": "Phyllodocida",
                        "family": "Nereididae", "genus": "Alitta"}}
    t = ms.build_table([d], lineage=lineage)
    assert list(t["accession"]) == ["MN000001.1", "MN000002.1"]
    r = t.iloc[0]
    assert r["length"] == 15720 and r["family"] == "Nereididae"
    assert r["seq_category"] == "NGS" and r["assembler_group"] == "Multi-assembly"
    assert r["assemblers"] == "MitoZ; NOVOPlasty" and r["journal"] == "Mitochondrial DNA B Resour"
    assert r["year"] == 2020
    r2 = t.iloc[1]
    assert r2["seq_category"] == "Unknown" and r2["assembler_group"] == "Not stated"  # no Assembly-Data block


def test_run_writes_table_and_report(tmp_path):
    d = tmp_path / "gb"
    d.mkdir()
    (d / "batch_0001.gb").write_text(REC)
    out = tmp_path / "out"
    ms.run([str(d)], str(out), lineage={"6346": {"phylum": "Annelida"}}, figures=False)
    t = pd.read_csv(out / "mitostats.csv")
    assert t.loc[0, "seq_category"] == "NGS"
    html = (out / "mitostats_report.html").read_text(encoding="utf-8")
    data = json.loads(html.split('<script id="g2t-data" type="application/json">')[1].split("</script>")[0])
    assert data["records"][0]["accession"] == "MN000001.1"
    assert Path(out / "counts_by_year.csv").exists()


def test_citations_fill_unpublished_journals(tmp_path):
    t = pd.DataFrame({"accession": ["A1.1", "A2.1", "A3.1"], "pub_status": ["unpublished", "published", "unpublished"],
                      "journal": ["", "Zool Stud", ""], "pub_year": [None, 2020, None], "year": [2019, 2020, 2018],
                      "year_source": ["submission", "publication", "submission"]})
    pd.DataFrame({"accession": ["A1.1", "A2.1"], "Ref1Journal": ["Mitochondrial DNA B Resour 5, 1-2 (2021)",
                                                                 "Other J 1, 1 (2022)"]}).to_csv(
        tmp_path / "cite_curate.csv", index=False)
    r = ms.apply_citations(t, tmp_path / "cite.csv")
    assert list(r["journal"]) == ["Mitochondrial DNA B Resour", "Zool Stud", ""]  # GenBank journal kept
    assert list(r["journal_source"]) == ["Crossref", "GenBank", ""]
    assert r.loc[0, "year"] == 2021 and r.loc[2, "pub_status"] == "unpublished"


@pytest.mark.parametrize(
    "a, b, same",
    [
        ("Mitochondrial DNA B Resour", "Mitochondrial DNA Part B", True),
        ("Mitochondrial DNA B Resour", "Mitochondrial DNA Part B: Resources", True),
        ("Mol Phylogenet Evol", "Molecular Phylogenetics and Evolution", True),
        ("Zool Stud", "Zoological Studies", True),
        ("Mitochondrial DNA A DNA Mapp Seq Anal", "Mitochondrial DNA Part B", False),
        ("Gene", "Genes", False),
        ("Diversity", "Diversity and Distributions", False),
        ("Diversity", "Mar Pollut Bull", False),
    ],
)
def test_same_journal(a, b, same):
    assert ms.same_journal(a, b) is same


def test_citation_journal_takes_genbank_spelling(tmp_path):
    t = pd.DataFrame({"accession": ["A0.1", "A0.2", "A1.1", "A2.1"],
                      "pub_status": ["published", "published", "published", "unpublished"],
                      "journal": ["Mitochondrial DNA", "Mitochondrial DNA", "Mitochondrial DNA B Resour", ""],
                      "pub_year": [2014, 2014, 2020, None], "year": [2014, 2014, 2020, 2019],
                      "year_source": ["publication"] * 3 + ["submission"]})
    pd.DataFrame({"accession": ["A2.1"], "Ref1Journal": ["Mitochondrial DNA Part B 6 (3), 1-2 (2021)"]}).to_csv(
        tmp_path / "c_curate.csv", index=False)
    # "Mitochondrial DNA" is commoner, but "Mitochondrial DNA B Resour" is the specific match
    assert ms.apply_citations(t, tmp_path / "c.csv").loc[3, "journal"] == "Mitochondrial DNA B Resour"


@pytest.mark.parametrize(
    "raw, group",
    [
        ("Burrows-Wheeler Aligner v. 0.7.17; bcftools v. 1.9", "Reference mapping"),  # human mtDNA style
        ("MIA v. 1.0; leeHom v. 1", "Reference mapping"),  # ancient DNA
        ("BWA v. 0.7; SPAdes v. 3.15", "SPAdes"),  # mapping then de novo: the assembler counts
        ("Gnenious v. 11", "Geneious"),
        ("genious prime v. 2022", "Geneious"),
        ("NOVAPlasty v. 4.2", "NOVOPlasty"),
        ("Mito-Z v. 2.4", "MitoZ"),
        ("CodonCode Aligner v. 9", "CodonCode"),
        ("MAFFT v. 7; MEGA v. 11", "Not stated"),  # alignment only
        ("x", "Not stated"),
        ("build 2025-03-24", "Not stated"),
    ],
)
def test_assembler_group_from_metazoa_spellings(raw, group):
    names, other = ms.assemblers(raw)
    assert ms.assembler_group(names, other, ms.mapping_tools(raw)) == group


def test_report_template_does_not_spread_record_arrays():
    """Math.min(...years) overflowed the call stack with 150,025 Metazoa records (report never drew)."""
    from importlib.resources import files

    tpl = files("g2t").joinpath("data/mitostats_report.html").read_text(encoding="utf-8")
    assert not re.search(r"Math\.(min|max)\(\s*\.\.\.", tpl)


def test_exclude_taxa():
    t = pd.DataFrame({"organism": ["Homo sapiens", "Homo sapiens neanderthalensis", "Homo naledi", "Pan troglodytes",
                                   "Alitta virens"],
                      "class": ["Mammalia"] * 4 + ["Polychaeta"], "genus": ["Homo"] * 3 + ["Pan", "Alitta"]})
    assert list(ms.exclude_taxa(t, ["Homo sapiens"])["organism"]) == ["Homo naledi", "Pan troglodytes", "Alitta virens"]
    assert list(ms.exclude_taxa(t, ["Mammalia"])["organism"]) == ["Alitta virens"]
    assert len(ms.exclude_taxa(t, ["Homo"])) == 2


def test_taxonomy_cache_is_shared_between_output_folders(tmp_path):
    d = tmp_path / "gb"
    d.mkdir()
    out = tmp_path / "old_report"
    out.mkdir()
    (out / "taxonomy_lineage.json").write_text('{"6346": {"phylum": "Annelida"}}')
    assert ms._lineage_cache([str(d)], out) == d / "taxonomy_lineage.json"
    assert (d / "taxonomy_lineage.json").exists()  # copied from the earlier report folder
    assert ms._lineage_cache([str(d)], tmp_path / "other_report") == d / "taxonomy_lineage.json"
