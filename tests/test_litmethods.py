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


def test_sra_runs_and_instrument_platform(tmp_path, monkeypatch):
    from g2t import litmethods as lm
    from g2t import mitostats as ms

    runs = ms.sra_runs(["BioProject:PRJEB1", "Sequence Read Archive:ERR1000001, ERR1000002,", "ERR1000003"], "reads SRR123456")
    assert runs == ["ERR1000001", "ERR1000002", "ERR1000003", "SRR123456"]
    tsv = ("run_accession\tinstrument_platform\tinstrument_model\tlibrary_strategy\n"
           "ERR1000001\tPACBIO_SMRT\tRevio\tWGS\nERR1000002\tILLUMINA\tIllumina NovaSeq 6000\tHi-C\n")
    monkeypatch.setattr(lm, "_get", lambda url: tsv.encode())
    got = lm.sra_methods(["ERR1000001", "ERR1000002"], lm.Cache(tmp_path / "c.json"))
    assert got["sra_seq_platforms"] == "PacBio" and got["sra_seq_category"] == "Long-read"  # Hi-C run ignored
    assert got["sra_complete"] == "yes"


def test_sra_partial_and_ambiguous(tmp_path, monkeypatch):
    from g2t import litmethods as lm

    tsv = ("run_accession\tinstrument_platform\tinstrument_model\tlibrary_strategy\n"
           "ERR1\tILLUMINA\tIllumina NovaSeq 6000\tWGS\nERR2\tILLUMINA\tIllumina NovaSeq 6000\tOTHER\n")
    monkeypatch.setattr(lm, "_get", lambda url: tsv.encode())
    c = lm.Cache(tmp_path / "c.json")
    part = lm.sra_methods(["ERR1", "ERR9"], c)  # ERR9 unknown to ENA -> partial, not complete
    assert part["sra_complete"] == "no" and part["sra_unresolved"] == "ERR9"
    assert lm.sra_methods(["ERR2"], c) == {}  # ambiguous strategy (Hi-C is often filed as OTHER): nothing inferred
    assert f"{lm.CACHE_VERSION}:r:ERR9" not in c.d  # unknown runs are asked again next time
    row = {"seq_category": "Long-read", "seq_platforms": "", "sra_runs": "ERR1", "assembler_group": "OATK"}
    assert not lm.needs_lookup(row)  # platform named in COMMENT is not overridden


def test_biosample_runs_via_ena(tmp_path, monkeypatch):
    from g2t import mitostats as ms2
    assert ms2.biosamples(["BioProject:PRJEB1", "BioSample:SAMEA112947850"]) == ["SAMEA112947850"]
    monkeypatch.setattr(lm, "_get", lambda url: b"run_accession\tinstrument_platform\tinstrument_model\tlibrary_strategy\n"
                        b"ERR1\tILLUMINA\tIllumina NovaSeq 6000\tWGS\n")
    cache = lm.Cache(tmp_path / "c.json")
    runs = lm.ena_sample_runs(["SAMEA112947850"], cache)
    assert runs == ["ERR1"]
    assert lm.sra_methods(runs, cache)["sra_seq_platforms"] == "Illumina"


def test_record_pubmed_ids_column():
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import Reference
    rec = SeqRecord(Seq("ACGT"), id="X1.1")
    ref = Reference()
    ref.pubmed_id = "123"
    rec.annotations["references"] = [ref]
    assert [str(r.pubmed_id) for r in rec.annotations["references"]] == ["123"]


def test_assembled_and_annotated_sentence():
    t = "The clean reads were assembled de novo using MitoZ v3.4 and the mitogenome was annotated with MITOS."
    assert lm.parse_methods(t)["lit_assemblers"] == "MitoZ"
    t2 = "The mitogenome was annotated with MITOS after assembly."
    assert not lm.parse_methods(t2).get("lit_assemblers")


# real wordings seen in the 2026-11 spot check of "cited" papers (ERGA-BGE genome papers)
_FOAM = ("Finally, the mitochondrial genome ( OZ239036.1 ) was assembled as a single circular contig of 16,836 bp "
         "using the FOAM pipeline v0.5 and included in the released assembly (GCA_965183895.1).")
_CAPTION = "The Hi-C contact map shows the 23 autosomes and the mitochondrial genome (GenBank accession: OZ199603.1 )."
_TREE = "Holothuria forskali ( NC_013884.1 ), Holothuria sanctori ( OZ199603.1 ), Bohadschia argus ( NC_061402.1 )"


def test_own_genome_paper_is_reports():
    assert lm.accession_role(_FOAM, "OZ239036.1")[0] == "reports"
    assert lm.accession_role(_CAPTION, "OZ199603.1")[0] == "reports"


def test_species_list_stays_cited():
    assert lm.accession_role(_TREE, "OZ199603.1")[0] == "cited"


def test_title_names_record_species():
    t = "ERGA-BGE reference genome of &lt;i&gt;Holothuria (Platyperona) sanctori&lt;/i&gt;"
    assert lm.title_names_organism(t, "Holothuria sanctori")
    assert not lm.title_names_organism("Complete mitochondrial genome of Holothuria atra", "Holothuria sanctori")
    assert lm.accession_role(_TREE, "OZ199603.1", t, "Holothuria sanctori")[0] == "reports"
    assert lm.accession_role(_TREE, "OZ199603.1", "Holothuria atra mitogenome", "Holothuria sanctori")[0] == "cited"


def test_accession_absent_from_text():
    assert lm.accession_role("A paper about nutrition.", "OZ410380.1")[0] == "absent"
