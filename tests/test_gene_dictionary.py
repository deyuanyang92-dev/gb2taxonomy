"""The gene dictionary: one packaged source, owned per run, and the gene lists of download and
reconcile kept in step with it."""

from __future__ import annotations

import logging
import re
import subprocess
import sys

import pandas as pd
import pytest

from g2t import ncbi_genes
from g2t.classify import (
    MatchConfig,
    default_gene_dict,
    load_path_dict,
    match_nuclear_genes,
    read_gene_dict_file,
)
from g2t.reconcile import GENE_FULL, voucher_core

# gene types that are not a downloadable locus of their own
_NOT_A_LOCUS = {"mtgenome", "plastome", "18-28s"}
# classify type -> ncbi_genes key, where the names differ
_DOWNLOAD_NAME = {"cox2": "COII", "cox3": "COIII", "cob": "CYTB", "ef-1": "EF1A", "its1-its2": "ITS"}
# downloadable, but deliberately not a gene type of their own
_DOWNLOAD_ONLY = {"5.8S": "inside its1-its2 / 18-28s records", "RPB2": "not classified yet"}


def test_default_gene_dict_is_a_fresh_copy():
    a = default_gene_dict()
    a["coi"].append("my spelling")
    a["new"] = ["x"]
    b = default_gene_dict()
    assert "my spelling" not in b["coi"] and "new" not in b


def test_configs_do_not_share_a_dictionary(tmp_path):
    f = tmp_path / "extra.txt"
    f.write_text("h3: my special histone marker\n", encoding="utf-8")
    custom = MatchConfig(gene_dict=load_path_dict(str(f)))
    plain = MatchConfig()
    text = "my special histone marker, partial sequence"
    assert "h3" in match_nuclear_genes(text, 400, custom)[0]
    assert "h3" not in match_nuclear_genes(text, 400, plain)[0]
    assert "my special histone marker" not in default_gene_dict()["h3"]


def test_load_path_dict_does_not_modify_base(tmp_path):
    f = tmp_path / "extra.csv"
    pd.DataFrame({"gene": ["coi", "newgene"], "synonyms": ["coxone,co-one", "new marker"]}).to_csv(f, index=False)
    base = default_gene_dict()
    before = {k: list(v) for k, v in base.items()}
    merged = load_path_dict(str(f), base)
    assert base == before
    assert {"coxone", "co-one"} <= set(merged["coi"]) and merged["newgene"] == ["new marker"]


@pytest.mark.parametrize(
    "name,content,error",
    [
        ("missing.txt", None, FileNotFoundError),
        ("empty.txt", "no colon here\n", ValueError),
        ("bad.csv", "a,b\n1,2\n", ValueError),
        ("dict.yaml", "coi: [x]\n", ValueError),
        ("nodict.py", "x = 1\n", ValueError),
    ],
)
def test_bad_user_dictionary_raises(tmp_path, name, content, error):
    """A user dictionary that cannot be used stops the run (it used to log and classify without it)."""
    f = tmp_path / name
    if content is not None:
        f.write_text(content, encoding="utf-8")
    with pytest.raises(error):
        read_gene_dict_file(str(f))


def _canon(name: str) -> str:
    return name.upper().replace("-", "").replace(".", "")


def test_every_gene_type_can_be_downloaded():
    downloadable = {_canon(g) for g in ncbi_genes.GENES}
    missing = [
        g for g in default_gene_dict() if g not in _NOT_A_LOCUS and _canon(_DOWNLOAD_NAME.get(g, g)) not in downloadable
    ]
    assert not missing, f"add these gene types to g2t/ncbi_genes.py GENES: {missing}"


def test_every_download_gene_is_a_gene_type():
    types = {_canon(_DOWNLOAD_NAME.get(g, g)) for g in default_gene_dict()}
    extra = [g for g in ncbi_genes.GENES if _canon(g) not in types and g not in _DOWNLOAD_ONLY]
    assert not extra, f"downloaded but never classified (add to data/gene_dict.yaml or _DOWNLOAD_ONLY): {extra}"


def test_every_single_locus_gene_type_is_a_voucher_gene_token():
    single = [g for g in default_gene_dict() if g not in _NOT_A_LOCUS and "-" not in g or g == "ef-1"]
    missing = [g for g in single if not GENE_FULL.fullmatch(g)]
    assert not missing, f"add to g2t/reconcile.py GENE_TOKENS: {missing}"


@pytest.mark.parametrize("gene", ["ND3", "nd4l", "NAD6", "ATP8"])
def test_all_mitochondrial_protein_genes_downloadable(gene):
    assert ncbi_genes.clause(gene).startswith("(")


def test_plastid_gene_after_voucher_is_stripped():
    # 'KUN 12345 matK' gave no voucher core before matK was a gene token
    assert voucher_core("KUN 12345 matK") == voucher_core("KUN12345") == "KUN12345"


def test_cli_logs_once_and_writes_classify_log(tmp_path):
    f = tmp_path / "final.csv"
    pd.DataFrame(
        [
            {
                "LocusID": "A",
                "ACCESSION": "A",
                "Length": 658,
                "MoleculeType": "DNA",
                "mol_type": "genomic DNA",
                "organelle": "mitochondrion",
                "Topology": "linear",
                "Definition": "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene; mitochondrial",
                "Organism": "Priapulus caudatus",
            }
        ]
    ).to_csv(f, index=False)
    out = tmp_path / "out"
    code = "import sys; from g2t._cli import main_classify; main_classify(sys.argv[1:])"
    p = subprocess.run(
        [sys.executable, "-c", code, "-i", str(f), "-o", str(out)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=True,
    )
    console = p.stdout + p.stderr
    lines = [ln for ln in console.splitlines() if "Round 1: Gene type matching" in ln]
    assert len(lines) == 1, console
    assert "Round 1: Gene type matching" in (out / "classify.log").read_text(encoding="utf-8")
    assert not (tmp_path / "extract_methods.log").exists()
    assert not logging.getLogger("g2t").handlers


def test_gene_dict_yaml_has_no_duplicate_types():
    from importlib.resources import files

    import yaml

    text = (files("g2t") / "data" / "gene_dict.yaml").read_text(encoding="utf-8")
    keys = re.findall(r"^  ([^\s:#][^:]*):\s*$", text, re.M)
    assert len(keys) == len(set(keys))
    assert set(yaml.safe_load(text)["gene_types"]) == set(default_gene_dict())


def test_compartment_lists_cover_the_dictionary_exactly():
    from g2t.genes import GENE_ORDER, MITO_TYPES, NUCLEAR_TYPES, PLASTID_TYPES

    groups = [set(MITO_TYPES), set(NUCLEAR_TYPES), set(PLASTID_TYPES)]
    assert sum(len(g) for g in groups) == len(GENE_ORDER) == len(set(GENE_ORDER)), "a type in two compartments"
    assert set(GENE_ORDER) == set(default_gene_dict())


def test_steps_use_the_shared_gene_lists():
    from g2t import curate, genes, organize, reconcile

    cfg = MatchConfig()
    assert cfg.mito_genes == set(genes.MITO_TYPES) and cfg.plastid_genes == set(genes.PLASTID_TYPES)
    assert organize.OrganizeConfig().gene_order == list(genes.GENE_ORDER)
    mito_cols = organize.OrganizeConfig().get_mito_genes()
    assert set(mito_cols) == set(genes.MITO_TYPES) - {"mtgenome"}
    assert curate.GENE_TOKENS == reconcile.GENE_TOKENS == genes.VOUCHER_GENE_TOKENS
    from g2t._classify.genedict import _DEFAULT_NUCLEAR_GENES

    assert set(genes.NUCLEAR_TYPES) == _DEFAULT_NUCLEAR_GENES


def test_classify_cli_defaults_are_the_shared_lists():
    from g2t import classify, genes

    seen = {}

    def fake_run(args):
        seen.update(vars(args))

    orig = classify._run_cli
    classify._run_cli = fake_run
    try:
        classify.main(["-i", "x.csv"])
    finally:
        classify._run_cli = orig
    assert seen["mtgenes_list"].split(",") == list(genes.MITO_TYPES)
    assert seen["ntgenes_list"].split(",") == list(genes.NUCLEAR_TYPES)
    assert seen["cpgenes_list"].split(",") == list(genes.PLASTID_TYPES)


def test_curate_strips_the_same_gene_names_as_reconcile():
    from g2t.curate import normalize_voucher

    assert normalize_voucher("ND4L_ZMMU_WS399") == "ZMMU_WS399"
    assert normalize_voucher("matK KUN 12345") == "KUN_12345"


def test_genome_records_fill_every_gene_column_of_their_compartment():
    from g2t import genes, organize

    inc = organize.OrganizeConfig().gene_includes
    assert inc["mtgenome"] == [g for g in genes.MITO_TYPES if g != "mtgenome"]
    assert inc["plastome"] == [g for g in genes.PLASTID_TYPES if g != "plastome"]
    assert set(inc["18-28s"]) <= set(genes.NUCLEAR_TYPES)


def test_explicit_mtgenome_includes_are_honoured_exactly():
    """Passing the six genes that used to be the stock list was silently widened to all 15."""
    from g2t.organize import OrganizeConfig, _mtgenome_children

    six = ["coi", "16s", "12s", "cob", "cox2", "cox3"]
    cfg = OrganizeConfig()
    cfg.gene_includes["mtgenome"] = list(six)
    assert _mtgenome_children(cfg) == six


def test_curate_searches_every_gene_column():
    from g2t import curate, genes

    assert list(genes.GENE_ORDER) == curate.DEFAULT_GENES


def test_download_dry_run_genes_are_downloadable_gene_types():
    for g in ncbi_genes.MITO_GENES + ncbi_genes.NUCLEAR_GENES + ncbi_genes.PLASTID_GENES:
        assert ncbi_genes.canon(g) == g


def test_length_ranges_name_only_known_gene_types():
    from g2t import genes

    assert set(MatchConfig().gene_length_ranges) <= set(genes.GENE_ORDER)


def test_histone3_spelling_is_h3(tmp_path):
    """'histone3 (H3)' (no space; 23 Enchytraeidae records, e.g. ON491619) was unmatched."""
    from g2t.classify import classify

    f = tmp_path / "final.csv"
    pd.DataFrame([{"LocusID": "ON491619", "ACCESSION": "ON491619", "Length": 328, "MoleculeType": "DNA",
                   "mol_type": "genomic DNA", "organelle": "", "Topology": "linear",
                   "Definition": "Hemienchytraeus wuhanensis voucher NH5 histone3 (H3) gene, partial cds",
                   "Organism": "Hemienchytraeus wuhanensis"}]).to_csv(f, index=False)
    classify(str(f), str(tmp_path / "out"))
    st = pd.read_csv(tmp_path / "out" / "record_status.csv")
    assert st.loc[0, "status"] == "assigned" and st.loc[0, "gene_type"] == "h3"
