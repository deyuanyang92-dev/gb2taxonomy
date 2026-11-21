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
