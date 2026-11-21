import json

from g2t import litmethods as lm

XML = ("<article><body>"
       "<sec><title>Introduction</title><p>Earlier work used PacBio and MitoZ.</p></sec>"
       "<sec><title>Materials and methods</title>"
       "<sec><title>Sequencing</title><p>The library was sequenced on Illumina NovaSeq. "
       "The mitochondrial genome was assembled using GetOrganelle and NOVOPlasty. "
       "Genes were annotated with MITOS.</p></sec></sec>"
       "<sec><title>Results</title><p>The mitogenome was submitted to GenBank (accession no. OP104125).</p></sec>"
       "</body><back><ref-list><ref>MitoZ: a toolkit</ref></ref-list></back></article>")


def test_sections_and_parse():
    meth, whole, found = lm.article_sections(XML)
    assert found and "PacBio" not in meth and "a toolkit" not in whole
    p = lm.parse_methods(meth)
    assert p["lit_seq_platforms"] == "Illumina"
    assert set(p["lit_assemblers"].split("; ")) == {"GetOrganelle", "NOVOPlasty"}


def test_role_reports_vs_cited():
    _, whole, _ = lm.article_sections(XML)
    assert lm.accession_role(whole, "OP104125.1")[0] == "reports"
    other = "Sequences were downloaded from GenBank (OP104125) as outgroup."
    assert lm.accession_role(other, "OP104125")[0] == "cited"


def test_exact_accession_token():
    rx = lm.accession_re("AB123456.1")
    assert rx.search("AB123456") and rx.search("(AB123456.2)")
    assert not rx.search("AB1234567") and not rx.search("XAB123456")


def test_prose_false_positives():
    p = lm.parse_methods("Reagents from Roche were used. Solid glass beads were added. A 454 bp band appeared.")
    assert p["lit_seq_platforms"] == ""


def test_search_follows_cursor(monkeypatch):
    pages = [{"resultList": {"result": [{"pmcid": "A"}]}, "nextCursorMark": "c1"},
             {"resultList": {"result": [{"pmcid": "B"}]}, "nextCursorMark": "c1"}]
    monkeypatch.setattr(lm, "_get", lambda url, retries=4: json.dumps(pages.pop(0)).encode())
    assert [h["pmcid"] for h in lm.search_accession("X1")] == ["A", "B"]


def test_damaged_cache_is_ignored(tmp_path):
    p = tmp_path / "c.json"
    p.write_text("{broken")
    c = lm.Cache(p)
    assert c.get("k", lambda: 1) == 1 and json.loads(p.read_text())
