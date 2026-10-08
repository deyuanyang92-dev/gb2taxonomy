"""Offline tests for g2t.download query building and g2t.ncbi_genes."""

import pytest

from g2t.download import DownloadOptions, build_query
from g2t.ncbi_genes import GENES, canon, clause


class TestGenes:
    @pytest.mark.parametrize("name,expected", [
        ("coi", "COI"), ("COX1", "COI"), ("co1", "COI"), ("18s", "18S"), ("SSU", "18S"),
        ("ef1alpha", "EF1A"), ("5.8S", "5.8S"), ("cob", "CYTB"), ("ITS2", "ITS"),
    ])
    def test_canon(self, name, expected):
        assert canon(name) == expected

    def test_unknown_gene(self):
        with pytest.raises(ValueError):
            canon("XYZ9")

    def test_clause_is_or_group(self):
        c = clause("COI")
        assert c.startswith("(") and c.endswith(")") and " OR " in c
        assert all(x in c for x in GENES["COI"])


class TestBuildQuery:
    def test_default_markers(self):
        q, tag = build_query("37891", DownloadOptions())
        assert tag == "markers"
        assert q.startswith("txid37891[Organism:exp]")
        for x in ("NOT wgs[filter]", "NOT biomol_mrna[PROP]", "NOT srcdb_refseq_model[PROP]", "NOT refseq[filter]"):
            assert x in q

    def test_all(self):
        q, tag = build_query("1", DownloadOptions(all=True))
        assert tag == "all" and q == "txid1[Organism:exp]"

    def test_mito_and_mitogenome(self):
        q, tag = build_query("1", DownloadOptions(mito=True))
        assert tag == "mito" and "mitochondrion[filter]" in q
        q, tag = build_query("1", DownloadOptions(mitogenome=True))
        assert tag == "mitogenome" and "10000:30000[SLEN]" in q

    def test_gene_and_length(self):
        q, tag = build_query("1", DownloadOptions(gene="coi,18s", minlen=500))
        assert tag == "gene-COI_18S-len500-max"
        assert "COX1[Gene]" in q and "18S[Title]" in q and "500:999999999[SLEN]" in q

    def test_includes_and_custom(self):
        opts = DownloadOptions(include_wgs=True, include_refseq=True, query="Russia[Country]", tag="x")
        q, tag = build_query("1", opts)
        assert "wgs[filter]" not in q and "NOT refseq[filter]" not in q
        assert "AND (Russia[Country])" in q and tag == "x"
