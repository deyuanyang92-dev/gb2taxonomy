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
            return [  # NCBI answers with the current accession.version, also for a bare accession
                {"AccessionVersion": a if "." in a else a + ".1", "Length": 15000 + i, "CreateDate": "2020/01/02",
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


def test_failed_summary_part_is_retried_next_run(tmp_path, fake, monkeypatch):
    fake["list"] = accs(6)
    real = dl._call
    calls = {"n": 0}

    def flaky(fn, parse, retries=5, **kw):
        if fn.__name__ == "esummary" and "MN000004.1" in kw["id"] and calls["n"] == 0:
            calls["n"] += 1
            raise RuntimeError("NCBI 500")
        return real(fn, parse, retries, **kw)

    monkeypatch.setattr(dl, "_call", flaky)
    opts = DownloadOptions(mitogenome=True, header_only=True, batch=3, workers=2)
    res = download("Annelida", str(tmp_path), opts, progress=False)
    with open(Path(res.out_dir) / "summary.tsv", encoding="utf-8") as fh:
        assert [r["accession"] for r in csv.DictReader(fh, delimiter="\t")] == accs(3)
    fake["summary_ids"].clear()
    download("Annelida", str(tmp_path), opts, progress=False)
    assert fake["summary_ids"] == accs(6)[3:]  # only the part that failed
    with open(Path(res.out_dir) / "summary.tsv", encoding="utf-8") as fh:
        assert [r["accession"] for r in csv.DictReader(fh, delimiter="\t")] == accs(6)


# ---------------------------------------------------------------- ideas taken from kingfisher


def test_summary_only_writes_table_and_no_records(tmp_path, fake):
    fake["list"] = accs(4)
    res = download("Annelida", str(tmp_path), DownloadOptions(mitogenome=True, summary_only=True), progress=False)
    out = Path(res.out_dir)
    assert out.name == "mitogenome-summary"
    assert not fake["header_calls"] and not fake["full_calls"]
    assert not list(out.glob("batch_*.gb")) and not (tmp_path / "_records.sqlite").exists()
    with open(out / "summary.tsv", encoding="utf-8") as fh:
        assert [r["accession"] for r in csv.DictReader(fh, delimiter="\t")] == accs(4)


def test_accession_list_download(tmp_path, fake):
    # unversioned accessions get the current version from esummary; given versions are kept
    res = download("", str(tmp_path), DownloadOptions(header_only=True),
                   accessions=["MN000001", "MN000002.1", "MN000002.1", ""], progress=False)
    assert res.count == 2 and res.fetched == 2
    assert Path(res.out_dir).name.startswith("list-") and Path(res.out_dir).name.endswith("-header")
    assert sorted(sum(fake["header_calls"], [])) == ["MN000001.1", "MN000002.1"]
    again = download("", str(tmp_path), DownloadOptions(header_only=True),
                     accessions=["MN000001", "MN000002.1"], progress=False)
    assert again.reused == 2 and again.fetched == 0  # rerun: nothing fetched


def test_read_accession_file(tmp_path):
    f = tmp_path / "a.csv"
    f.write_text('accession,organism\n# comment\nMN000001.1,Alitta\n"OQ123456",x\nPP_1\nNC_012345.1\n')
    assert dl.read_accession_file(str(f)) == ["MN000001.1", "OQ123456", "NC_012345.1"]


def test_progress_line_has_rate_and_eta():
    line = dl._progress(10, 100, 500, 5000, 120.0)
    assert line.startswith("[10/100] 5,000 records  2,500/min  ETA 18 min")
    assert dl._progress(100, 100, 500, 50000, 600.0).split("  ")[2] == "done"


def test_more_workers_with_api_key(monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))  # no saved key
    assert dl.default_workers() == 3 and dl.default_workers("KEY") == 10
    monkeypatch.setenv("NCBI_API_KEY", "x")
    assert dl.default_workers() == 10


def test_saved_api_key_is_used(monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    assert dl.saved_api_key() is None
    (tmp_path / "g2t").mkdir()
    (tmp_path / "g2t" / "credentials").write_text("email=\napi_key = abc123\n")
    assert dl.saved_api_key() == "abc123" and dl.default_workers() == 10
