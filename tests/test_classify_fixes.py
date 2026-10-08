"""Fixes from the code review of classify.py / utils.py (L-03, L-06, L-07, L-10, L-13, L-14, L-15, L-16,
L-17, C1-C4)."""

import os
import subprocess
import sys

import pandas as pd
import pytest

from g2t.classify import MatchConfig, classify, clean_definition, clean_length, process_row, smart_match
from g2t.utils import length_in_range, parse_interval


def _row(acc, length, definition, organelle="", **kw):
    r = {"LocusID": acc, "ACCESSION": acc, "Version": f"{acc}.1", "Length": length, "MoleculeType": "DNA",
         "mol_type": "genomic DNA", "organelle": organelle, "Topology": "linear", "Definition": definition,
         "Organism": "Priapulus caudatus", "organism": "Priapulus caudatus"}
    r.update(kw)
    return r


def _write(tmp_path, rows, name="final.csv"):
    f = tmp_path / name
    pd.DataFrame(rows).to_csv(f, index=False)
    return f


@pytest.mark.parametrize("raw,val", [("550 bp", 550), (550.0, 550), ("550.0", 550), ("1,234 bp", 1234),
                                     ("", 0), (None, 0)])
def test_l03_clean_length(raw, val):
    assert clean_length(raw) == val


def test_l03_missing_length_row_does_not_inflate_others(tmp_path):
    rows = [_row("A", "550 bp", "Priapulus caudatus 18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S "
                 "ribosomal RNA gene, internal transcribed spacer 2, and 28S ribosomal RNA gene, partial sequence"),
            _row("B", "", "Priapulus caudatus something")]
    out = tmp_path / "o"
    classify(str(_write(tmp_path, rows)), str(out), MatchConfig())
    st = pd.read_csv(out / "record_status.csv", dtype=str).fillna("")
    assert st.loc[st.LocusID == "A", "gene_type"].iloc[0] != "18-28s"


@pytest.mark.parametrize("definition,organelle,bad", [
    ("Homo sapiens 16S rRNA methyltransferase (rsmA) gene, complete cds", "mitochondrion", "16s"),
    ("histone H3 lysine 4 demethylase gene, partial cds", "", "h3"),
    ("NADH dehydrogenase subunit 1 (ND1) gene and its flanking sequence, partial cds", "", "its1-its2"),
    ("Priapulus caudatus elongation factor-1 beta gene, partial cds", "", "ef-1"),
    ("Priapulus caudatus isolate CO2 16S ribosomal RNA gene, partial sequence; mitochondrial",
     "mitochondrion", "cox2"),
])
def test_l06_false_positives(definition, organelle, bad):
    res = process_row(_row("X", 900, definition, organelle), MatchConfig())
    assert bad not in res["gene_type"].split(",")


def test_l06_true_positives_kept():
    cfg = MatchConfig()
    assert process_row(_row("A", 330, "Priapulus caudatus histone H3 gene, partial cds, and its upstream region"),
                       cfg)["gene_type"] == "h3"
    assert "its1-its2" in process_row(_row("B", 700, "Priapulus caudatus internal transcribed spacer 1, "
                                                     "partial sequence (ITS1)"), cfg)["gene_type"]
    assert process_row(_row("C", 658, "Priapulus caudatus voucher WS399 cytochrome c oxidase subunit I (COI) "
                                      "gene, partial cds; mitochondrial", "mitochondrion"), cfg)["gene_type"] == "coi"


def test_l06_clean_definition():
    assert "co2" not in clean_definition("Priapulus caudatus isolate CO2 16S ribosomal RNA gene")
    assert "its" not in clean_definition("gene and its flanking region").split()
    assert "its1" in clean_definition("internal transcribed spacer 1 (ITS1)")


def test_l07_deterministic():
    d = ("Priapulus caudatus mitochondrial small subunit ribosomal RNA and large subunit ribosomal RNA genes, "
         "partial sequence; mitochondrial")
    code = ("from g2t.classify import process_row, MatchConfig;"
            f"print(process_row({{'LocusID':'X','Length':900,'Definition':{d!r},'organelle':'mitochondrion',"
            "'MoleculeType':'DNA','Topology':'linear'}, MatchConfig())['gene_type'])")
    outs = {subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                           env={**os.environ, "PYTHONHASHSEED": str(seed)}).stdout.strip() for seed in range(6)}
    assert len(outs) == 1, outs


def test_l15_synonym_starting_with_punctuation():
    assert smart_match("(coi)", "cytochrome oxidase (coi) gene")
    assert smart_match("-cox1", "x -cox1 y")


def test_l16_parse_interval():
    assert parse_interval("none") == (None, None) and parse_interval("all") == (None, None)
    assert parse_interval("150:") == (150, None) and parse_interval(":500") == (None, 500)
    with pytest.raises(ValueError):
        parse_interval("500:150")
    assert length_in_range(100, "all")


def test_l17_mtgenome_lower_bound_from_config():
    d = "Priapulus caudatus mitochondrion, complete genome"
    cfg = MatchConfig()
    cfg.gene_length_ranges["mtgenome"] = "2000:500000"
    assert process_row(_row("A", 2500, d, "mitochondrion"), cfg)["gene_type"] == "mtgenome"


def test_l10_c2_stale_outputs_removed(tmp_path):
    out = tmp_path / "o"
    rows = [_row("A", 658, "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene; mitochondrial",
                 "mitochondrion"),
            _row("B", 700, "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene; mitochondrial",
                 "mitochondrion")]
    classify(str(_write(tmp_path, rows)), str(out), MatchConfig())
    classify(str(_write(tmp_path, rows[:1], "f2.csv")), str(out), MatchConfig(), if_recheck=False)
    all_df = pd.read_csv(out / "assigned_genes_types_all.csv", dtype=str)
    assert list(all_df["LocusID"]) == ["A"]
    assert not (out / "assigned_genes_types2.csv").exists()


def test_l10_all_file_written_without_recheck(tmp_path):
    out = tmp_path / "o"
    rows = [_row("A", 658, "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene; mitochondrial",
                 "mitochondrion")]
    r = classify(str(_write(tmp_path, rows)), str(out), MatchConfig(), if_recheck=False)
    assert r.success and r.rows == 1 and (out / "assigned_genes_types_all.csv").exists()


def test_c3_all_records_filtered(tmp_path):
    out = tmp_path / "o"
    rows = [_row("A", "50 bp", "Priapulus caudatus COI gene"), _row("B", "60 bp", "Priapulus caudatus 18S")]
    r = classify(str(_write(tmp_path, rows)), str(out), MatchConfig())
    assert r.success and r.rows == 0
    st = pd.read_csv(out / "record_status.csv", dtype=str)
    assert set(st["status"]) == {"filtered"}


def test_c4_thousands_separator(tmp_path):
    out = tmp_path / "o"
    rows = [_row("A", "1,234 bp", "Priapulus caudatus 18S ribosomal RNA gene, partial sequence")]
    classify(str(_write(tmp_path, rows)), str(out), MatchConfig())
    st = pd.read_csv(out / "record_status.csv", dtype=str)
    assert st.loc[0, "status"] == "assigned"


def test_c1_cli_writes_record_status(tmp_path):
    rows = [_row("A", 658, "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene; mitochondrial",
                 "mitochondrion"), _row("B", 100, "short thing")]
    f = _write(tmp_path, rows)
    out = tmp_path / "cli"
    subprocess.run([sys.executable, "-m", "g2t.classify", "-i", str(f), "-o", str(out)], check=True,
                   capture_output=True, cwd=tmp_path)
    st = pd.read_csv(out / "record_status.csv", dtype=str)
    assert set(st["status"]) == {"assigned", "filtered"}


def test_l14_column_names_with_spaces_kept(tmp_path):
    out = tmp_path / "o"
    rows = [_row("A", 658, "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene; mitochondrial",
                 "mitochondrion", **{"Assembly Method": "SPAdes"})]
    classify(str(_write(tmp_path, rows)), str(out), MatchConfig())
    assert "Assembly Method" in pd.read_csv(out / "assigned_genes_types_all.csv").columns
