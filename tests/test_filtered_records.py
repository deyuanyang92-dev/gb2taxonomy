"""Every input record must be accounted for: assigned, unmatched or filtered (with a reason)."""

import pandas as pd

from g2t.classify import MatchConfig, classify


def _row(acc, length, definition, **kw):
    r = {"LocusID": acc, "ACCESSION": acc, "Version": f"{acc}.1", "Length": f"{length} bp",
         "MoleculeType": "DNA", "mol_type": "genomic DNA", "organelle": "", "Definition": definition,
         "Organism": "Priapulus caudatus", "organism": "Priapulus caudatus"}
    r.update(kw)
    return r


def _input(tmp_path):
    rows = [
        _row("A1", 658, "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene, partial cds; mitochondrial"),
        _row("A2", 1800, "Priapulus caudatus 18S ribosomal RNA gene, partial sequence"),
        _row("A3", 120, "Priapulus caudatus FOXA3 factor gene, partial cds"),            # too short
        _row("A4", 600, "Priapulus caudatus FOXA3 factor gene, partial cds"),            # unmatched
        _row("A5", 90000, "Priapulus caudatus chromosome X, partial"),                   # too long
        _row("A6", "", "Priapulus caudatus something"),                                  # no length
        _row("A2", 1800, "Priapulus caudatus 18S ribosomal RNA gene, partial sequence"),  # duplicate LocusID
    ]
    f = tmp_path / "final.csv"
    pd.DataFrame(rows).to_csv(f, index=False)
    return f


def test_filtered_records_written_with_reasons(tmp_path):
    out = tmp_path / "out"
    r = classify(str(_input(tmp_path)), str(out), MatchConfig())
    assert r.success
    filt = pd.read_csv(out / "filtered_records.csv", dtype=str).fillna("")
    reasons = dict(zip(filt["ACCESSION"], filt["reason"]))
    assert "A3" in reasons and "length 120 bp outside 150:50000" in reasons["A3"]
    assert "A5" in reasons and "outside 150:50000" in reasons["A5"]
    assert "A6" in reasons and "length missing" in reasons["A6"]
    assert (filt["reason"].str.contains("duplicate LocusID")).sum() == 1
    assert set(filt["step"]) <= {"classify round 1", "classify round 2"}


def test_record_status_accounts_for_every_input_row(tmp_path):
    out = tmp_path / "out"
    classify(str(_input(tmp_path)), str(out), MatchConfig())
    st = pd.read_csv(out / "record_status.csv", dtype=str).fillna("")
    assert len(st) == 7                                   # one row per input row, duplicates included
    status = dict(zip(st["input_row"], st["status"]))
    assert status["0"] == "assigned" and status["1"] == "assigned"
    assert status["2"] == "filtered" and status["4"] == "filtered" and status["5"] == "filtered"
    assert status["3"] == "unmatched"
    assert status["6"] == "filtered"                      # duplicate of A2
    assert (st.loc[st["status"] != "assigned", "reason"] != "").all()
    assert st.loc[st["input_row"] == "0", "gene_type"].iloc[0] == "coi"


def test_wide_length_range_moves_short_records_to_unmatched(tmp_path):
    out = tmp_path / "out"
    classify(str(_input(tmp_path)), str(out), MatchConfig(global_length_range="1:1000000"))
    st = pd.read_csv(out / "record_status.csv", dtype=str).fillna("")
    status = dict(zip(st["input_row"], st["status"]))
    assert status["2"] == "unmatched" and status["4"] == "unmatched"
    assert status["5"] == "filtered"                      # still no length


def test_rerun_overwrites_previous_lists(tmp_path):
    out = tmp_path / "out"
    classify(str(_input(tmp_path)), str(out), MatchConfig())
    classify(str(_input(tmp_path)), str(out), MatchConfig())
    filt = pd.read_csv(out / "filtered_records.csv", dtype=str)
    assert len(filt) == 4                                 # A3, A5, A6, duplicate A2 - not doubled
