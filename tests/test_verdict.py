import pandas as pd

from g2t.verdict import combine


def _df(**kw):
    return pd.DataFrame([kw])


def test_two_layers_far_is_contaminant():
    v = combine(tree=_df(accession="A.1", gene="coi", organism="x", flag="outside_family"),
                blast=_df(accession="A.1", gene="coi", conflict_rank="family", gate="actionable", lca_rank="family"))
    assert v.verdict.tolist() == ["likely_contaminant"]


def test_two_layers_near_is_misidentified():
    v = combine(kmer=_df(accession="A.1", gene="16s", organism="x", flag="near_other_genus"),
                blast=_df(accession="A.1", gene="16s", conflict_rank="genus", gate="actionable", lca_rank="genus"))
    assert v.verdict.tolist() == ["likely_misidentified"]


def test_single_layer_suspect_and_numt_and_clean():
    v = combine(tree=_df(accession="A.1", gene="coi", organism="x", flag="outside_family"))
    assert v.verdict.tolist() == ["suspect"]
    v = combine(trans=_df(accession="B.1", gene="coi", flag="numt_suspect"))
    assert v.verdict.tolist() == ["numt_suspect"]
    v = combine(blast=_df(accession="C.1", gene="coi", conflict_rank="", gate="actionable", lca_rank="genus",
                          agrees=True))
    assert v.verdict.tolist() == ["clean"]
    v = combine(blast=_df(accession="D.1", gene="coi", conflict_rank="", gate="no_close_reference", lca_rank=""))
    assert v.verdict.tolist() == ["unchecked"]


def test_genes_are_not_pooled():
    v = combine(kmer=_df(accession="A.1", gene="coi", organism="x", flag="remote"),
                tree=_df(accession="A.1", gene="16s", organism="x", flag="inside_other_genus"))
    assert sorted(zip(v.gene, v.verdict)) == [("16s", "suspect"), ("coi", "suspect")]


def test_kmer_and_tree_are_one_family():
    v = combine(kmer=_df(accession="A.1", gene="coi", organism="x", flag="remote"),
                tree=_df(accession="A.1", gene="coi", organism="x", flag="outside_family"))
    assert v.verdict.tolist() == ["suspect"]
    v = combine(kmer=_df(accession="A.1", gene="coi", organism="x", flag="remote"),
                tree=_df(accession="A.1", gene="coi", organism="x", flag="outside_family"), strict=False)
    assert v.verdict.tolist() == ["likely_contaminant"]


def test_far_tree_vs_near_blast_disagree():
    v = combine(tree=_df(accession="A.1", gene="coi", organism="x", flag="outside_family"),
                blast=_df(accession="A.1", gene="coi", conflict_rank="species", gate="actionable"))
    assert v.verdict.tolist() == ["likely_misidentified"]
    assert "layers_disagree" in v.reasons[0]


def test_geneless_evidence_applies_to_all_markers():
    v = combine(kmer=_df(accession="A.1", gene="", organism="x", flag="near_other_genus"),
                blast=pd.DataFrame([dict(accession="A.1", gene=g, conflict_rank="genus", gate="actionable")
                                    for g in ("coi", "16s")]))
    assert v.verdict.tolist() == ["likely_misidentified"] * 2


def test_family_agreement_is_not_clean():
    v = combine(blast=_df(accession="E.1", gene="16s", conflict_rank="", gate="actionable", agrees=True,
                          resolved=False, confirmed_rank="family"))
    assert v.verdict.tolist() == ["coarse_agree"] and v.confirmed_rank.tolist() == ["family"]
    v = combine(blast=_df(accession="E.1", gene="16s", conflict_rank="", gate="actionable", agrees=True,
                          resolved=True, confirmed_rank="species"))
    assert v.verdict.tolist() == ["clean"]


def test_destinations_must_agree_for_likely_misidentified():
    v = combine(kmer=_df(accession="A.1", gene="coi", organism="x", flag="near_other_genus", nearest_genus="Beta"),
                blast=_df(accession="A.1", gene="coi", conflict_rank="genus", gate="actionable", lca_rank="genus",
                          lca_taxon="Gamma"))
    assert v.verdict.tolist() == ["suspect"] and "destinations_disagree" in v.reasons[0]
    v = combine(kmer=_df(accession="A.1", gene="coi", organism="x", flag="near_other_genus", nearest_genus="Beta"),
                blast=_df(accession="A.1", gene="coi", conflict_rank="species", gate="actionable", lca_rank="species",
                          lca_taxon="Beta alba"))
    assert v.verdict.tolist() == ["likely_misidentified"]


def test_weak_blast_conflict_alone_is_suspect_and_never_escalates():
    v = combine(blast=_df(accession="A.1", gene="coi", conflict_rank="genus", gate="actionable", weak=True))
    assert v.verdict.tolist() == ["suspect"]
    v = combine(kmer=_df(accession="A.1", gene="coi", organism="x", flag="near_other_genus", nearest_genus="Beta"),
                blast=_df(accession="A.1", gene="coi", conflict_rank="genus", gate="actionable", weak=True))
    assert v.verdict.tolist() == ["suspect"]


def test_no_hit_row_is_unchecked_not_missing():
    v = combine(blast=_df(accession="A.1", gene="coi", gate="no_hit", note="no BLAST hit"))
    assert v.verdict.tolist() == ["unchecked"]


def test_reverse_complement_is_suspect_not_clean():
    v = combine(blast=_df(accession="A.1", gene="coi", organism="x", gate="actionable", agrees=True,
                          resolved=True, confirmed_rank="species", reverse=True))
    assert v.verdict.tolist() == ["suspect"]
    assert "reverse_complement" in v.reasons[0]


def _inv(rows):
    return pd.DataFrame(rows, columns=["accession", "gene", "organism", "family", "taxonomy"])


def test_every_record_gets_a_row_without_evidence():
    inv = _inv([("A1.1", "COI", "Nereis sp.", "Nereididae", "Eukaryota; Metazoa; Annelida")])
    v = combine(inventory=inv, taxon="Nereididae")
    assert len(v) == 1 and v.iloc[0]["flag"] == "NA" and "BLAST" in v.iloc[0]["why"]


def test_non_metazoa_and_other_family_are_check():
    inv = _inv([("B1.1", "COI", "Aeromonas sp.", "", "Bacteria; Proteobacteria"),
                ("B2.1", "COI", "Syllis sp.", "Syllidae", "Eukaryota; Metazoa; Annelida")])
    v = combine(inventory=inv, taxon="Nereididae").set_index("accession")
    assert v.loc["B1", "flag"] == "CHECK" and "不是后生动物" in v.loc["B1", "why"]
    assert v.loc["B2", "flag"] == "CHECK" and "研究类群" in v.loc["B2", "why"]


def test_ncbi_unverified_and_blast_no_hit_are_check():
    inv = _inv([("C1.1", "COI", "X y", "Nereididae", "Eukaryota; Metazoa"),
                ("C2.1", "COI", "X z", "Nereididae", "Eukaryota; Metazoa")])
    kmer = pd.DataFrame([dict(accession="C1", gene="COI", flag="ncbi_unverified", organism="X y")])
    blast = pd.DataFrame([dict(accession="C2", gene="COI", gate="no_hit", conflict_rank="")])
    v = combine(kmer=kmer, blast=blast, inventory=inv, taxon="Nereididae", ran=("blast",)).set_index("accession")
    assert v.loc["C1", "flag"] == "CHECK" and "ncbi_unverified" in v.loc["C1", "why"]
    assert v.loc["C2", "flag"] == "CHECK" and "无命中" in v.loc["C2", "why"]
