"""--header-only: records without sequence (efetch seq_start=1&seq_stop=1) plus lengths from esummary."""

import csv
import re
from pathlib import Path

import pytest

import g2t.download as dl
from g2t.download import DownloadOptions, build_query, download


def _rec(acc, length):
    return (
        f"LOCUS       {acc.split('.')[0]}                  {length} bp    DNA     linear   INV 01-JAN-2020\n"
        f"DEFINITION  record {acc}.\nACCESSION   {acc.split('.')[0]}\nVERSION     {acc}\n"
        "ORIGIN\n        1 a\n//\n"
    )


def accs(n):
    return [f"MN{i:06d}.1" for i in range(1, n + 1)]


@pytest.fixture
def fake(monkeypatch):
    state = {"list": [], "header_calls": [], "full_calls": [], "summary_ids": []}

    def _call(fn, parse, retries=5, **kw):
        name = fn.__name__
        if kw.get("db") == "taxonomy":
            return (
                {"IdList": ["6340"]}
                if name == "esearch"
                else [{"ScientificName": "Annelida", "Rank": "phylum", "Lineage": "x"}]
            )
        if name == "esearch":
            return {"Count": str(len(state["list"])), "WebEnv": "w", "QueryKey": "1"}
        if name == "esummary":
            ids = kw["id"].split(",")
            state["summary_ids"] += ids
            return [
                {"AccessionVersion": a, "Length": 15000 + i, "CreateDate": "2020/01/02",
                 "UpdateDate": "2021/03/04", "TaxId": 6341, "Title": f"Worm {a} mitochondrion, complete genome"}
                for i, a in enumerate(ids)
            ]
        if kw.get("rettype") == "acc":
            return "\n".join(state["list"][kw["retstart"] : kw["retstart"] + kw["retmax"]]) + "\n"
        raise AssertionError(name)

    monkeypatch.setattr(dl, "_call", _call)

    def _efetch_text(ids, header_only=False):
        (state["header_calls"] if header_only else state["full_calls"]).append(list(ids))
        return "".join(_rec(a, 1 if header_only else 15000) for a in ids)

    monkeypatch.setattr(dl, "_efetch_text", _efetch_text)
    monkeypatch.setattr(dl.time, "sleep", lambda s: None)
    return state


def test_header_only_changes_the_tag():
    _, tag = build_query("6340", DownloadOptions(mitogenome=True, header_only=True))
    assert tag == "mitogenome-header"
    q_full, _ = build_query("6340", DownloadOptions(mitogenome=True))
    q_head, _ = build_query("6340", DownloadOptions(mitogenome=True, header_only=True))
    assert q_full == q_head  # same records, only the content differs


def test_header_only_fetches_headers_and_writes_lengths(tmp_path, fake):
    fake["list"] = accs(7)
    res = download("Annelida", str(tmp_path), DownloadOptions(mitogenome=True, header_only=True, batch=3),
                   progress=False)
    assert res.fetched == 7 and not fake["full_calls"]
    assert sum(len(c) for c in fake["header_calls"]) == 7
    out = Path(res.out_dir)
    assert out.name == "mitogenome-header"
    with open(out / "summary.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["accession"] for r in rows] == accs(7)
    assert rows[0]["length"] == "15000" and rows[0]["create_date"] == "2020-01-02"
    assert rows[0]["taxid"] == "6341" and rows[0]["update_date"] == "2021-03-04"
    # headers live in their own store: a later full download must not take them for full records
    assert (tmp_path / "_records_header.sqlite").exists()
    assert not (tmp_path / "_records.sqlite").exists()


def test_full_download_ignores_header_batches(tmp_path, fake):
    fake["list"] = accs(4)
    download("Annelida", str(tmp_path), DownloadOptions(mitogenome=True, header_only=True), progress=False)
    res = download("Annelida", str(tmp_path), DownloadOptions(mitogenome=True), progress=False)
    assert res.fetched == 4 and res.reused == 0
    text = (Path(res.out_dir) / "batch_0001.gb").read_text()
    assert re.findall(r"(\d+) bp", text) == ["15000"] * 4


def test_rerun_queries_only_new_records(tmp_path, fake):
    fake["list"] = accs(5)
    download("Annelida", str(tmp_path), DownloadOptions(mitogenome=True, header_only=True), progress=False)
    fake["header_calls"].clear()
    fake["summary_ids"].clear()
    fake["list"] = accs(7)  # two new mitogenomes at NCBI
    res = download("Annelida", str(tmp_path), DownloadOptions(mitogenome=True, header_only=True), progress=False)
    assert res.reused == 5 and res.fetched == 2
    assert sum(fake["header_calls"], []) == accs(7)[5:]
    assert fake["summary_ids"] == accs(7)[5:]
    with open(Path(res.out_dir) / "summary.tsv", encoding="utf-8") as fh:
        assert [r["accession"] for r in csv.DictReader(fh, delimiter="\t")] == accs(7)
