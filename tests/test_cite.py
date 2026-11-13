"""g2t-cite: Crossref matching of GenBank records still marked Unpublished (no network: fake search)."""

import json

import pandas as pd
import pytest

from g2t import cite


@pytest.mark.parametrize(
    "authors, names",
    [
        ("Radashevsky,V.I., Pankova,V.V. and Lee,B.", ["radashevsky", "pankova", "lee"]),
        ("Smith,A.", ["smith"]),
        ("Muller,J. and Mueller,K.", ["muller", "mueller"]),
        ("", []),
    ],
)
def test_surnames(authors, names):
    assert cite.surnames(authors) == names


def test_title_sim_ignores_case_markup_and_punctuation():
    a = "The complete mitochondrial genome of Alitta virens (Annelida: Nereididae)"
    b = "The complete mitochondrial genome of <i>Alitta virens</i> (Annelida, Nereididae)"
    assert cite.title_sim(a, b) == 1.0
    assert cite.title_sim(a, "Phylogeny of Sabellida") < 0.5


@pytest.mark.parametrize(
    "ts, ov, tier",
    [(0.95, 1.0, "high"), (0.8, 0.5, "medium"), (0.7, 1.0, "check"), (0.95, 0.0, "check"), (0.6, 0.2, "low"),
     (0.55, 1.0, "low")],  # changed title without the taxon: too weak
)
def test_grade(ts, ov, tier):
    assert cite.grade(ts, ov) == tier


def test_one_shared_common_surname_is_not_medium():
    assert cite.grade(0.77, 0.5, shared=1) == "low"  # "four mitogenomes of Oligochaeta" vs a bird paper
    assert cite.grade(0.77, 0.5, shared=2) == "medium"
    assert cite.grade(0.55, 1.0, shared=3, taxon_named=True) == "check"


@pytest.mark.parametrize(
    "journals, unpub",
    [
        (["Unpublished", "Submitted (01-JAN-2020) X"], True),
        (["Zool Stud 65, e37 (2026)", "Submitted (01-JAN-2020) X"], False),
        (["Submitted (01-JAN-2020) X"], True),
        ([], False),
    ],
)
def test_is_unpublished(journals, unpub):
    assert cite.is_unpublished(journals) is unpub


ITEM = {
    "DOI": "10.1080/23802359.2020.1775147",
    "type": "journal-article",
    "title": ["The complete mitochondrial genome of <i>Alitta virens</i>"],
    "author": [{"family": "Smith", "given": "Anna B."}, {"family": "Lee", "given": "Bo"}],
    "short-container-title": ["Mitochondrial DNA B Resour"],
    "volume": "5", "issue": "3", "page": "2370-2371",
    "published-print": {"date-parts": [[2020, 7]]},
}
FIGURE = {**ITEM, "type": "component", "DOI": "10.6084/fig1", "title": ["Figure 1 from: The complete mitochondrial"
                                                                         " genome of Alitta virens"]}


def test_crossref_record_in_genbank_style():
    assert cite.genbank_journal(ITEM) == "Mitochondrial DNA B Resour 5 (3), 2370-2371 (2020)"
    assert cite.cr_authors(ITEM) == "Smith,A.B. and Lee,B."


def test_best_match_skips_figure_components():
    m = cite.best_match("The complete mitochondrial genome of Alitta virens", "Smith,A. and Lee,B.", 2020,
                        lambda t, a, y: [FIGURE, ITEM])
    assert m["tier"] == "high" and m["doi"] == ITEM["DOI"]
    assert cite.best_match("x", "Smith,A.", None, lambda t, a, y: [FIGURE]) == {"tier": "none"}


REC = """LOCUS       MN000001                   1 bp    DNA     linear   INV 01-JAN-2021
DEFINITION  Alitta virens mitochondrion, complete genome.
ACCESSION   MN000001
VERSION     MN000001.1
SOURCE      mitochondrion Alitta virens
  ORGANISM  Alitta virens
            Eukaryota; Metazoa; Annelida.
REFERENCE   1  (bases 1 to 1)
  AUTHORS   Smith,A. and Lee,B.
  TITLE     The complete mitochondrial genome of Alitta virens
  JOURNAL   Unpublished
REFERENCE   2  (bases 1 to 1)
  AUTHORS   Smith,A.
  TITLE     Direct Submission
  JOURNAL   Submitted (01-JAN-2020) Somewhere
FEATURES             Location/Qualifiers
     source          1
                     /organism="Alitta virens"
ORIGIN
        1 a
//
"""


def test_run_groups_records_caches_and_writes_curate_sheet(tmp_path):
    d = tmp_path / "gb"
    d.mkdir()
    published = REC.replace("MN000001", "MN000003").replace("JOURNAL   Unpublished", "JOURNAL   Zool Stud 1, 1 (2020)")
    (d / "b.gb").write_text(REC + REC.replace("MN000001", "MN000002") + published)
    calls = []

    def search(t, a, y):
        calls.append((t, a, y))
        return [ITEM]

    out = tmp_path / "res" / "cite.xlsx"
    cand = cite.run([str(d)], str(out), search=search, delay=0, progress=False)
    assert calls == [("The complete mitochondrial genome of Alitta virens", "smith", 2020)]  # one query per paper
    assert len(cand) == 1 and cand.loc[0, "n_records"] == 2 and cand.loc[0, "tier"] == "high"
    cur = pd.read_excel(out, sheet_name="Curate")
    assert list(cur.columns) == cite.CURATE_COLUMNS
    assert list(cur["accession"]) == ["MN000001.1", "MN000002.1"]
    assert cur.loc[0, "Ref1Journal"] == "Mitochondrial DNA B Resour 5 (3), 2370-2371 (2020)"
    assert json.loads((out.parent / "crossref_cache.json").read_text())  # rerun: no new query
    cite.run([str(d)], str(out), search=search, delay=0, progress=False)
    assert len(calls) == 1


def test_table_input_with_ref_columns(tmp_path):
    f = tmp_path / "final.csv"
    pd.DataFrame({
        "accession": ["OQ1.1", "OQ2.1"], "organism": ["A b", "C d"],
        "Ref1Journal": ["Unpublished", "Zool Stud 1, 1 (2020)"], "Ref1Title": ["Paper one", "Paper two"],
        "Ref1Authors": ["Wu,X.", "Li,Y."], "Ref2Journal": ["Submitted (02-FEB-2019) X", ""],
        "Ref2Title": ["Direct Submission", ""], "Ref2Authors": ["Wu,X.", ""],
    }).to_csv(f, index=False)
    g = cite.unpublished_groups(cite.read_records([str(f)]))
    assert list(g["accessions"]) == ["OQ1.1"] and g.loc[0, "submit_year"] == 2019


def test_check_rows_are_not_in_curate_sheet(tmp_path):
    other = {**ITEM, "title": ["Mitochondrial genome of Alitta virens and the phylogeny of Nereididae"]}
    cand = cite.match_groups(
        pd.DataFrame([{"title": "The complete mitochondrial genome of Alitta virens", "authors": "Smith,A. and Lee,B.",
                       "submit_year": 2020, "n_records": 1, "accessions": "MN1.1", "organisms": "Alitta virens"}]),
        lambda t, a, y: [other], delay=0, progress=False)
    assert cand.loc[0, "tier"] == "check"  # same authors, changed title
    assert cite.curate_table(cand).empty


@pytest.mark.parametrize(
    "gb, cr, orgs, ok",
    [
        ("The complete mitochondrion genome of Megascolecidae sp.", "The Complete Mitochondrion Genome of Ecklonia "
         "maxima", ["Megascolecidae sp. XX"], False),
        ("Mitogenome of Amynthas heteropodus", "The complete mitochondrial genome of Amynthas deogyusanensis",
         ["Amynthas heteropodus"], False),
        ("Mitogenome of Amynthas righii", "Complete mitochondrial genome of the earthworm Amynthas righii",
         ["Amynthas righii"], True),
        ("Mitogenome of Hediste sinesimplex", "Description of the mitogenome of H. sinesimplex",
         ["Hediste sinesimplex"], True),
        ("Mitogenomes of three earthworms", "Mitogenomes of three earthworms", ["Amynthas righii"], None),
    ],
)
def test_taxon_must_agree(gb, cr, orgs, ok):
    assert cite.taxon_agrees(gb, cr, cite.taxon_names(orgs)) is ok


def test_crossref_title_cleaned():
    assert cite.clean_title("A new pheretimoid earthworm of the genus <i>Metaphire</i> Sims &amp;amp; Easton,\n"
                            "    1972") == "A new pheretimoid earthworm of the genus Metaphire Sims & Easton, 1972"
    assert cite.genbank_journal({"short-container-title": ["ZK"], "container-title": ["ZooKeys"],
                                 "volume": "1231", "page": "293-310", "issued": {"date-parts": [[2025]]}}) == \
        "ZooKeys 1231, 293-310 (2025)"
