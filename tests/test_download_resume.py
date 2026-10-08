"""Offline tests of g2t.download resume logic (Entrez calls are faked)."""

import json
import re
from pathlib import Path

import pytest

import g2t.download as dl
from g2t.download import DownloadOptions, download


class FakeNCBI:
    def __init__(self):
        self.db: dict[str, list[str]] = {}   # query substring -> accessions
        self.default: list[str] = []

    def result_for(self, term):
        for k, v in self.db.items():
            if k in term:
                return v
        return self.default


def _rec(acc):
    return (f"LOCUS       {acc.split('.')[0]}                  60 bp    DNA     linear   INV 01-JAN-2020\n"
            f"DEFINITION  record {acc}.\nACCESSION   {acc.split('.')[0]}\nVERSION     {acc}\n"
            "ORIGIN\n        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n//\n")


@pytest.fixture
def fake(monkeypatch):
    f = FakeNCBI()
    state = {}

    def _call(fn, parse, retries=5, **kw):
        name = fn.__name__
        if kw.get("db") == "taxonomy":
            return {"IdList": ["37891"]} if name == "esearch" else [
                {"ScientificName": "Priapulidae", "Rank": "family", "Lineage": "x"}]
        if name == "esearch":
            accs = f.result_for(kw["term"])
            if kw.get("usehistory"):
                state["last"] = accs
                return {"Count": str(len(accs)), "WebEnv": "w", "QueryKey": "1"}
            return {"Count": str(len(accs))}
        if kw.get("rettype") == "acc":
            return "\n".join(state["last"][kw["retstart"]: kw["retstart"] + kw["retmax"]]) + "\n"
        return "".join(_rec(a) for a in kw["id"].split(","))

    monkeypatch.setattr(dl, "_call", _call)
    monkeypatch.setattr(dl.time, "sleep", lambda s: None)
    return f


def accs(prefix, n):
    return [f"{prefix}{i:06d}.1" for i in range(1, n + 1)]


def ids_in(d):
    out = []
    for f in sorted(Path(d).glob("batch_*.gb")):
        out += re.findall(r"^VERSION\s+(\S+)", f.read_text(), re.M)
    return out


def test_growing_database_redownloads_shifted_batches(tmp_path, fake):
    fake.default = accs("AB", 400)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    fake.default = accs("ZZ", 3) + accs("AB", 400)          # 3 new records at the head of the list
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    got = ids_in(tmp_path / "markers")
    assert sorted(got) == sorted(fake.default) and res.failed == 0


def test_unchanged_list_resumes_without_downloading(tmp_path, fake):
    fake.default = accs("AB", 400)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    res = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    assert res.skipped == 2 and res.downloaded == 0


def test_different_queries_get_different_directories(tmp_path, fake):
    fake.db = {"Russia": accs("RU", 400), "Norway": accs("NO", 450)}
    r1 = download("Priapulidae", str(tmp_path), DownloadOptions(query="Russia[Country]"), progress=False)
    r2 = download("Priapulidae", str(tmp_path), DownloadOptions(query="Norway[Country]"), progress=False)
    assert r1.tag != r2.tag
    assert all(i.startswith("NO") for i in ids_in(r2.out_dir)) and len(ids_in(r2.out_dir)) == 450


def test_include_flags_change_the_tag(tmp_path, fake):
    fake.db = {"NOT wgs[filter]": accs("AB", 400)}
    fake.default = accs("WG", 400)
    r1 = download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    r2 = download("Priapulidae", str(tmp_path), DownloadOptions(include_wgs=True), progress=False)
    assert r1.tag == "markers" and r2.tag == "incl-wgs"
    assert all(i.startswith("WG") for i in ids_in(r2.out_dir))


def test_same_tag_different_query_is_refused(tmp_path, fake):
    fake.db = {"Russia": accs("RU", 10), "Norway": accs("NO", 10)}
    download("Priapulidae", str(tmp_path), DownloadOptions(query="Russia[Country]", tag="x"), progress=False)
    with pytest.raises(ValueError, match="different query"):
        download("Priapulidae", str(tmp_path), DownloadOptions(query="Norway[Country]", tag="x"), progress=False)


def test_stale_batches_removed_when_batch_size_grows(tmp_path, fake):
    fake.default = accs("AB", 400)
    download("Priapulidae", str(tmp_path), DownloadOptions(batch=100), progress=False)
    download("Priapulidae", str(tmp_path), DownloadOptions(batch=200), progress=False)
    got = ids_in(tmp_path / "markers")
    assert len(got) == len(set(got)) == 400
    assert sorted(p.name for p in (tmp_path / "markers").glob("batch_*.gb")) == ["batch_0001.gb", "batch_0002.gb"]


def test_dry_run_writes_nothing(tmp_path, fake):
    fake.default = accs("AB", 400)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    mf = tmp_path / "markers" / "manifest.json"
    before = mf.read_text()
    fake.default = accs("AB", 405)
    download("Priapulidae", str(tmp_path), DownloadOptions(), dry_run=True)
    assert mf.read_text() == before
    assert not (tmp_path / "mito").exists()
    download("Priapulidae", str(tmp_path), DownloadOptions(mito=True), dry_run=True)
    assert not (tmp_path / "mito").exists()


def test_manifest_records_completion(tmp_path, fake):
    fake.default = accs("AB", 5)
    download("Priapulidae", str(tmp_path), DownloadOptions(), progress=False)
    m = json.loads((tmp_path / "markers" / "manifest.json").read_text())
    assert m["complete"] is True and m["count"] == 5 and len(m["accessions_sha1"]) == 40
