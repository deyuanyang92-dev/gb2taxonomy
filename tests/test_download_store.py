"""Offline tests of the accession-level record store used by g2t.download (Entrez calls are faked)."""

import json
import re
from pathlib import Path

import pytest

import g2t.download as dl
from g2t.download import DownloadOptions, build_query, download
from g2t.recstore import RecordStore, iter_records


def _rec(acc, length=60):
    return (f"LOCUS       {acc.split('.')[0]}                  {length} bp    DNA     linear   INV 01-JAN-2020\n"
            f"DEFINITION  record {acc}.\nACCESSION   {acc.split('.')[0]}\nVERSION     {acc}\n"
            "ORIGIN\n        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n//\n")


def accs(prefix, n, ver=1):
    return [f"{prefix}{i:06d}.{ver}" for i in range(1, n + 1)]


@pytest.fixture
def fake(monkeypatch):
    state = {"list": [], "fetched": [], "fail_once": set(), "mdat": [], "terms": [], "cur": []}

    def _call(fn, parse, retries=5, **kw):
        name = fn.__name__
        if kw.get("db") == "taxonomy":
            return {"IdList": ["37891"]} if name == "esearch" else [
                {"ScientificName": "Priapulidae", "Rank": "family", "Lineage": "x"}]
        if name == "esearch":
            lst = state["mdat"] if "[MDAT]" in kw["term"] else state["list"]
            state["terms"].append(kw["term"])
            state["cur"] = lst
            return {"Count": str(len(lst)), "WebEnv": "w", "QueryKey": "1"}
        if kw.get("rettype") == "acc":
            return "\n".join(state["cur"][kw["retstart"]: kw["retstart"] + kw["retmax"]]) + "\n"
        ids = kw["id"].split(",")
        state["fetched"] += ids
        return "".join(_rec(a) for a in ids)

    monkeypatch.setattr(dl, "_call", _call)

    def _efetch_text(ids):
        def efetch():
            pass
        return _call(efetch, False, db="nuccore", id=",".join(ids), rettype="gbwithparts", retmode="text")

    monkeypatch.setattr(dl, "_efetch_text", _efetch_text)
    monkeypatch.setattr(dl.time, "sleep", lambda s: None)
    return state


def ids_in(d):
    out = []
    for f in sorted(Path(d).glob("batch_*.gb")):
        out += re.findall(r"^VERSION\s+(\S+)", f.read_text(), re.M)
    return out


def test_iter_records_splits_text():
    txt = _rec("AB000001.1") + _rec("AB000002.3")
    got = list(iter_records(txt.splitlines(keepends=True)))
    assert [a for a, _ in got] == ["AB000001.1", "AB000002.3"]
    assert got[1][1].endswith("//\n")


def test_store_roundtrip(tmp_path):
    s = RecordStore(tmp_path / "r.sqlite")
    s.put_many([("AB000001.1", _rec("AB000001.1"))])
    assert "AB000001.1" in s and "AB000001.2" not in s
    assert s.get("AB000001.1") == _rec("AB000001.1")
    assert s.missing(["AB000001.1", "AB000002.1"]) == ["AB000002.1"]


def test_growing_database_fetches_only_new_records(tmp_path, fake):
    fake["list"] = accs("AB", 400)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    fake["fetched"].clear()
    fake["list"] = accs("ZZ", 3) + accs("AB", 400)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert sorted(fake["fetched"]) == accs("ZZ", 3)
    assert res.fetched == 3 and res.reused == 400 and res.failed == 0
    assert sorted(ids_in(tmp_path / "markers")) == sorted(fake["list"])


def test_other_selection_reuses_records(tmp_path, fake):
    fake["list"] = accs("AB", 300)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    fake["fetched"].clear()
    fake["list"] = accs("AB", 300)[:120]
    res = download("Priapulidae", str(tmp_path), DownloadOptions(mito=True), progress=False)
    assert fake["fetched"] == [] and res.reused == 120
    assert len(ids_in(tmp_path / "mito")) == 120


def test_shared_store_across_taxa(tmp_path, fake):
    store = tmp_path / "cache" / "records.sqlite"
    fake["list"] = accs("AB", 50)
    download("Priapulidae", str(tmp_path / "t1"), DownloadOptions(), progress=False, store=str(store))
    fake["fetched"].clear()
    download("Priapulidae", str(tmp_path / "t2"), DownloadOptions(), progress=False, store=str(store))
    assert fake["fetched"] == []


def test_new_version_is_fetched_and_reported(tmp_path, fake):
    fake["list"] = accs("AB", 5)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    fake["fetched"].clear()
    fake["list"] = accs("AB", 4) + ["AB000005.2", "CD000001.1"]
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert sorted(fake["fetched"]) == ["AB000005.2", "CD000001.1"]
    assert (res.new, res.updated, res.removed) == (1, 1, 0)
    ch = (tmp_path / "markers" / "changes.tsv").read_text()
    assert "updated\tAB000005.2\tAB000005.1" in ch and "new\tCD000001.1" in ch
    assert "AB000005.1" not in ids_in(tmp_path / "markers")


def test_removed_records_reported(tmp_path, fake):
    fake["list"] = accs("AB", 5)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    fake["list"] = accs("AB", 3)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert res.removed == 2 and len(ids_in(tmp_path / "markers")) == 3


def test_legacy_batch_files_are_imported(tmp_path, fake):
    old = tmp_path / "markers"
    old.mkdir()
    (old / "batch_0001.gb").write_text("".join(_rec(a) for a in accs("AB", 10)))
    fake["list"] = accs("AB", 12)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert sorted(fake["fetched"]) == accs("AB", 12)[10:]
    assert res.reused == 10 and res.fetched == 2


def test_dry_run_reports_cache_without_writing(tmp_path, fake):
    fake["list"] = accs("AB", 10)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    mf = (tmp_path / "markers" / "manifest.json").read_text()
    fake["list"] = accs("AB", 15)
    fake["fetched"].clear()
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), dry_run=True)
    assert fake["fetched"] == [] and res.reused == 10 and res.to_fetch == 5
    assert (tmp_path / "markers" / "manifest.json").read_text() == mf


def test_markers_exclude_huge_records_by_default():
    q, tag = build_query("1", DownloadOptions())
    assert "AND (1:100000[SLEN] OR (mitochondrion[filter]" in q and tag == "markers"
    q, _ = build_query("1", DownloadOptions(include_large=True))
    assert "[SLEN]" not in q
    q, _ = build_query("1", DownloadOptions(maxlen=5000))
    assert "1:5000[SLEN]" in q and "1:100000[SLEN]" not in q
    q, _ = build_query("1", DownloadOptions(all=True))
    assert "[SLEN]" not in q


def test_parallel_workers_give_same_result(tmp_path, fake):
    fake["list"] = accs("AB", 1000)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(batch=50, workers=4), progress=False)
    got = ids_in(tmp_path / "markers")
    assert got == fake["list"] and res.failed == 0


def test_query_change_is_flagged(tmp_path, fake):
    fake["list"] = accs("AB", 5)
    download("Priapulidae", str(tmp_path), DownloadOptions(include_large=True, tag="x"), progress=False)
    m = tmp_path / "x" / "manifest.json"
    m.write_text(m.read_text().replace("NOT refseq[filter]", "NOT refseq[filter] AND old"))
    fake["list"] = accs("AB", 3)
    with pytest.raises(ValueError, match="different query"):
        download("Priapulidae", str(tmp_path), DownloadOptions(include_large=True, tag="x"), progress=False)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert res.query_changed is False
    d = tmp_path / "markers" / "manifest.json"
    d.write_text(d.read_text().replace("NOT refseq[filter]", "NOT refseq[filter] AND old"))
    fake["list"] = accs("AB", 2)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert res.query_changed and res.removed == 1
    assert "removed(query changed)" in (tmp_path / "markers" / "changes.tsv").read_text()


def test_mitogenome_query_has_no_upper_length_and_takes_assembly_titles():
    q, tag = build_query("1", DownloadOptions(mitogenome=True))
    assert tag == "mitogenome" and "10000:999999999[SLEN]" in q and "30000" not in q
    assert '"genome assembly"[Title]' in q and "100000[SLEN]" not in q


def test_failing_record_is_isolated_by_splitting(tmp_path, fake, monkeypatch):
    bad = "AB000007.1"
    orig = dl._efetch_text

    def _efetch_text(ids):
        if bad in ids:
            raise OSError("IncompleteRead")
        return orig(ids)

    monkeypatch.setattr(dl, "_efetch_text", _efetch_text)
    fake["list"] = accs("AB", 20)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(batch=10), progress=False)
    got = ids_in(tmp_path / "markers")
    assert bad not in got and len(got) == 19 and res.failed == 1
    m = (tmp_path / "markers" / "manifest.json").read_text()
    assert '"complete": false' in m


def test_efetch_text_decompresses_gzip(monkeypatch):
    import gzip as gz
    import io

    class Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    seen = {}

    def urlopen(req, timeout):
        seen["enc"] = req.headers.get("Accept-encoding")
        seen["timeout"] = timeout
        return Resp(gz.compress(_rec("AB000001.1").encode()))

    monkeypatch.setattr(dl.urllib.request, "urlopen", urlopen)
    assert dl._efetch_text(["AB000001.1"]) == _rec("AB000001.1")
    assert seen == {"enc": "gzip", "timeout": dl.TIMEOUT}


def test_unchanged_rerun_keeps_batch_files(tmp_path, fake):
    fake["list"] = accs("AB", 30)
    download("Priapulidae", str(tmp_path), DownloadOptions(batch=10), progress=False)
    f = tmp_path / "markers" / "batch_0001.gb"
    f.write_text(f.read_text() + "")          # same content
    mt = f.stat().st_mtime_ns
    download("Priapulidae", str(tmp_path), DownloadOptions(batch=10), progress=False)
    assert f.stat().st_mtime_ns == mt
    fake["list"] = accs("AB", 31)
    download("Priapulidae", str(tmp_path), DownloadOptions(batch=10), progress=False)
    assert len(ids_in(tmp_path / "markers")) == 31


def test_download_skips_report_queries_unless_asked(tmp_path, fake, monkeypatch):
    fake["list"] = accs("AB", 5)
    terms = []
    orig = dl.count
    monkeypatch.setattr(dl, "count", lambda t: terms.append(t) or orig(t))
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert len(terms) == 1                                   # only the selection itself
    terms.clear()
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False, report=True)
    assert len(terms) > 10 and res.per_gene and res.composition
    terms.clear()
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), dry_run=True)
    assert len(terms) > 10


def _set_date(out, iso):
    m = out / "manifest.json"
    d = json.loads(m.read_text())
    d["date"] = iso
    m.write_text(json.dumps(d))


def test_since_auto_fetches_only_modified_records(tmp_path, fake):
    fake["list"] = accs("AB", 20)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    _set_date(tmp_path / "markers", "2026-01-10")
    fake["fetched"].clear()
    fake["terms"].clear()
    fake["mdat"] = ["AB000003.2", "CD000001.1"]          # one updated, one new
    res = download("Priapulidae", str(tmp_path), DownloadOptions(since="auto"), progress=False)
    assert any("2026/01/07:3000[MDAT]" in t for t in fake["terms"])          # last date minus 3 days
    assert sorted(fake["fetched"]) == ["AB000003.2", "CD000001.1"]
    got = ids_in(tmp_path / "markers")
    assert "AB000003.1" not in got and "AB000003.2" in got and "CD000001.1" in got and len(got) == 21
    assert (res.new, res.updated, res.removed, res.since) == (1, 1, 0, "2026/01/07")
    m = json.loads((tmp_path / "markers" / "manifest.json").read_text())
    assert m["since"] == "2026/01/07" and m["last_full_check"] is not None


def test_since_without_earlier_download_falls_back_to_full_list(tmp_path, fake):
    fake["list"] = accs("AB", 5)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(since="auto"), progress=False)
    assert res.since == "" and res.fetched == 5


def test_no_store_mode_is_incremental_without_sqlite(tmp_path, fake):
    fake["list"] = accs("AB", 12)
    download("Priapulidae", str(tmp_path), DownloadOptions(no_store=True, batch=5), progress=False)
    assert not list(tmp_path.glob("*.sqlite"))
    assert sorted(ids_in(tmp_path / "markers")) == accs("AB", 12)
    fake["fetched"].clear()
    fake["list"] = accs("AB", 11) + ["AB000012.2", "CD000001.1"]
    res = download("Priapulidae", str(tmp_path), DownloadOptions(no_store=True, batch=5), progress=False)
    assert sorted(fake["fetched"]) == ["AB000012.2", "CD000001.1"]
    got = ids_in(tmp_path / "markers")
    assert sorted(got) == sorted(fake["list"]) and len(got) == len(set(got))
    assert (res.reused, res.updated, res.new) == (11, 1, 1)
    fake["list"] = accs("AB", 3)                          # withdrawn records are removed from the files
    download("Priapulidae", str(tmp_path), DownloadOptions(no_store=True, batch=5), progress=False)
    assert sorted(ids_in(tmp_path / "markers")) == accs("AB", 3)
    assert not list(tmp_path.glob("*.sqlite"))


def test_no_store_with_since(tmp_path, fake):
    fake["list"] = accs("AB", 8)
    download("Priapulidae", str(tmp_path), DownloadOptions(no_store=True), progress=False)
    _set_date(tmp_path / "markers", "2026-02-01")
    fake["fetched"].clear()
    fake["mdat"] = ["EF000001.1"]
    download("Priapulidae", str(tmp_path), DownloadOptions(no_store=True, since="2026-03-01"), progress=False)
    assert fake["fetched"] == ["EF000001.1"]
    assert any("2026/03/01:3000[MDAT]" in t for t in fake["terms"])
    assert len(ids_in(tmp_path / "markers")) == 9


def test_cli_validates_arguments(capsys):
    from g2t.download import main
    for argv in (["-t", "X", "--since", "yesterday"], ["-t", "X", "--store", "s.db", "--no-store"],
                 ["-t", "X", "-w", "0"], ["-t", "X", "--minlen", "-5"]):
        with pytest.raises(SystemExit) as e:
            main(argv)
        assert e.value.code == 2
    assert "auto" in capsys.readouterr().err


def test_cli_quiet_and_since_date(tmp_path, fake, capsys):
    from g2t.download import main
    fake["list"] = accs("AB", 3)
    assert main(["-t", "X", "-o", str(tmp_path), "-q"]) == 0
    assert "[1/1]" not in capsys.readouterr().err
    fake["mdat"] = []
    assert main(["-t", "X", "-o", str(tmp_path), "--since", "2026/01/02", "-q"]) == 0
    assert any("2026/01/02:3000[MDAT]" in t for t in fake["terms"])
