"""Regression tests for extract.py fixes L-01, L-02, L-04, L-05, L-19, L-20, L-21, L-24.

All inputs are tiny hand-written GenBank records; no test touches the network
(Bio.Entrez is monkeypatched and ete3 is masked in every test).
"""

from __future__ import annotations

import importlib
import importlib.util
import io
import json
import logging
import os
import re
import sys

import pandas as pd
import pytest
from Bio import Entrez, SeqIO

# import_module: ``g2t.extract`` may be shadowed by the re-exported function in ``g2t/__init__``
ex = importlib.import_module("g2t.extract")
LOGGER_NAME = "g2t.extract"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _origin(n: int) -> str:
    seq = ("acgt" * (n // 4 + 1))[:n]
    lines = []
    for i in range(0, n, 60):
        chunk = seq[i:i + 60]
        blocks = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {blocks}")
    return "\n".join(lines)


def _qual_lines(quals: dict[str, object]) -> str:
    out = []
    for k, v in quals.items():
        vals = v if isinstance(v, list) else [v]
        for item in vals:
            out.append(f'                     /{k}="{item}"')
    return "\n".join(out)


def make_record(
    acc: str,
    *,
    n: int = 60,
    organism: str = "Priapulus caudatus",
    taxon: str | None = "6232",
    db_xrefs: list[str] | None = None,
    quals: dict[str, object] | None = None,
    source: bool = True,
    extra_sources: list[dict[str, object]] | None = None,
    assembly: dict[str, str] | None = None,
    pubmed: str | None = None,
    length_token: str | None = None,
) -> str:
    """Return one syntactically valid GenBank record as text."""
    xrefs = db_xrefs if db_xrefs is not None else ([f"taxon:{taxon}"] if taxon else [])
    q: dict[str, object] = {"organism": organism, "mol_type": "genomic DNA"}
    if xrefs:
        q["db_xref"] = xrefs
    q.update(quals or {})
    lines = [
        f"LOCUS       {acc:<16} {(length_token or str(n)):>11} bp    DNA     linear   INV 12-AUG-2025",
        f"DEFINITION  {organism} cytochrome oxidase subunit I (COI) gene, partial cds.",
        f"ACCESSION   {acc}",
        f"VERSION     {acc}.1",
        "KEYWORDS    .",
        f"SOURCE      {organism}",
        f"  ORGANISM  {organism}",
        "            Eukaryota; Metazoa; Priapulida; Priapulidae.",
        f"REFERENCE   1  (bases 1 to {n})",
        "  AUTHORS   Doe,J.",
        "  TITLE     Direct Submission",
        "  JOURNAL   Submitted (07-AUG-2025) Somewhere",
    ]
    if pubmed:
        lines.append(f"  PUBMED    {pubmed}")
    if assembly:
        lines.append("COMMENT     ##Assembly-Data-START##")
        for k, v in assembly.items():
            lines.append(f"            {k} :: {v}")
        lines.append("            ##Assembly-Data-END##")
    lines.append("FEATURES             Location/Qualifiers")
    if source:
        lines.append(f"     source          1..{n}")
        lines.append(_qual_lines(q))
        for extra in extra_sources or []:
            lines.append("     source          1..10")
            lines.append(_qual_lines(extra))
    lines.append(f"     gene            1..{n}")
    lines.append('                     /gene="COI"')
    lines.append("ORIGIN")
    lines.append(_origin(n))
    lines.append("//")
    return "\n".join(lines) + "\n"


def write(path, *records: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(records), encoding="utf-8")
    return str(path)


def read_str_csv(path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def parse_one(text: str):
    return SeqIO.read(io.StringIO(text), "genbank")


# ---------------------------------------------------------------------------
# fixtures: no network, no ete3, fresh taxonomy state
# ---------------------------------------------------------------------------
class FakeEntrez:
    """Stand-in for Bio.Entrez.efetch/read (db='taxonomy')."""

    def __init__(self, lineages: dict[str, list[tuple]] | None = None, fail: bool = False):
        # lineages: taxid -> [(rank, name), ...] (ancestors only), plus ("self", (rank, name))
        self.lineages = lineages or {}
        self.fail = fail
        self.calls: list[list[str]] = []
        self.seen_email: list[object] = []
        self.seen_api_key: list[object] = []

    def efetch(self, db, id=None, **kw):
        assert db == "taxonomy"
        if self.fail:
            raise OSError("simulated network failure")
        ids = id.split(",") if isinstance(id, str) else [str(i) for i in id]
        self.calls.append(ids)
        self.seen_email.append(Entrez.email)
        self.seen_api_key.append(Entrez.api_key)
        return _Handle(ids)

    def read(self, handle, **kw):
        recs = []
        for tid in handle.ids:
            if tid not in self.lineages:
                continue  # NCBI silently omits unknown ids
            ancestors, own = self.lineages[tid]
            recs.append({
                "TaxId": tid,
                "ScientificName": own[1],
                "Rank": own[0],
                "AkaTaxIds": [],
                "LineageEx": [
                    {"TaxId": str(i), "ScientificName": name, "Rank": rank}
                    for i, (rank, name) in enumerate(ancestors, 1)
                ],
            })
        return recs


class _Handle:
    def __init__(self, ids):
        self.ids = ids

    def close(self):
        pass


def _blocked(*args, **kwargs):
    raise OSError("network access is blocked in tests")


@pytest.fixture(autouse=True)
def isolated_taxonomy(monkeypatch):
    monkeypatch.setattr(ex, "_taxonomy", ex.TaxonomyService())
    # Reproduce Python 3.13: ete3 is installed (spec found) but `import ete3` fails (no `cgi`).
    monkeypatch.setitem(sys.modules, "ete3", None)  # `import ete3` -> ModuleNotFoundError
    real_find_spec = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec",
                        lambda name, *a, **k: object() if name == "ete3" else real_find_spec(name, *a, **k))
    monkeypatch.setattr(Entrez, "efetch", _blocked)
    monkeypatch.setattr(Entrez, "email", None)
    monkeypatch.setattr(Entrez, "api_key", None)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    # Make a failing network attempt fast even if the implementation retries.
    monkeypatch.setattr(Entrez, "sleep_between_tries", 0)


def install_fake_entrez(monkeypatch, fake: FakeEntrez) -> FakeEntrez:
    monkeypatch.setattr(Entrez, "efetch", fake.efetch)
    monkeypatch.setattr(Entrez, "read", fake.read)
    return fake


# ---------------------------------------------------------------------------
# L-01  --batch parallel mode
# ---------------------------------------------------------------------------
class TestL01Batch:
    @pytest.mark.parametrize("stream", [True, False])
    def test_batch_workers_run_and_final_csv_matches_non_batch(self, tmp_path, stream):
        f1 = write(
            tmp_path / "a" / "x.gb",
            make_record("TST000001", quals={"isolate": "I1"}, assembly={"Assembly Name": "asm1"}),
            make_record("TST000002", quals={"specimen_voucher": "V2"}),
        )
        # same basename in another directory -> must not overwrite the first
        f2 = write(
            tmp_path / "b" / "x.gb",
            make_record("TST000003", quals={"strain": "S3", "lat_lon": "1 N 2 E"}, pubmed="12345"),
        )
        out_b = tmp_path / "batch"
        summary = ex.process_batch_mode(
            [f1, f2], str(out_b), max_tasks=2, stream_mode=stream, include_taxonomy=False)
        statuses = [(os.path.basename(os.path.dirname(r.file_path)), r.status) for r in summary.file_results]
        assert all(s == ex.FileStatus.SUCCESS for _, s in statuses), \
            [(r.file_path, r.status, r.error_message) for r in summary.file_results]
        assert summary.total_metadata_extracted == 3

        out_n = tmp_path / "plain"
        out_n.mkdir()
        ex.process_all_files([f1, f2], str(out_n), include_taxonomy=False, stream_mode=stream)

        final_b = read_str_csv(out_b / "final.csv")
        final_n = read_str_csv(out_n / "final.csv")
        assert list(final_b["ACCESSION"]) == ["TST000001", "TST000002", "TST000003"]
        pd.testing.assert_frame_equal(final_b, final_n)  # identical columns, order and values

        # same-basename inputs were kept in separate per-file directories
        subdirs = [d for d in os.listdir(out_b) if d.endswith("_metadata")]
        assert len(subdirs) == 2

    def test_extract_batch_true_produces_final_csv(self, tmp_path):
        d = tmp_path / "in"
        write(d / "one.gb", make_record("TST000011"))
        write(d / "two.gb", make_record("TST000012"), make_record("TST000013"))
        out = tmp_path / "out"
        res = ex.extract([str(d)], str(out), stream=True, batch=True, max_tasks=2)
        assert res.success
        assert res.output_file == str(out / "final.csv")
        assert res.rows == 3
        assert sorted(read_str_csv(out / "final.csv")["ACCESSION"]) == ["TST000011", "TST000012", "TST000013"]


# ---------------------------------------------------------------------------
# L-02  mid-file parse failure, reconciliation, pipeline report
# ---------------------------------------------------------------------------
class TestL02Failures:
    def _broken_file(self, tmp_path):
        return write(
            tmp_path / "trunc.gb",
            make_record("TST000021"),
            make_record("TST000022", length_token="abc"),
            make_record("TST000023"),
        )

    def test_midfile_parse_error_is_not_success(self, tmp_path, caplog):
        fp = self._broken_file(tmp_path)
        out = tmp_path / "out"
        out.mkdir()
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            summary = ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=False)
        fr = summary.file_results[0]
        assert fr.status in (ex.FileStatus.PARTIAL, ex.FileStatus.FAILED)
        assert fr.status != ex.FileStatus.SUCCESS
        assert fr.metadata_extracted == 1
        assert fr.status == ex.FileStatus.PARTIAL  # one record was recovered
        assert fr.error_message
        assert summary.files_succeeded == 0 and summary.files_partial == 1
        # reconciliation against the '//' count is logged and recorded in the report
        text = " ".join(r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING)
        assert "trunc.gb" in text and re.search(r"\b3\b.*\b1\b", text), text
        assert any(e.phase == "parse" for e in fr.record_errors)

    def test_truncated_download_is_flagged_by_terminator_count(self, tmp_path, caplog):
        """Biopython happily parses a cut-off last record; the '//' count exposes it."""
        text = make_record("TST000024", n=200) + make_record("TST000025", n=200)[:-120]
        fp = write(tmp_path / "cut.gb", text)
        out = tmp_path / "out"
        out.mkdir()
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            summary = ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=False)
        fr = summary.file_results[0]
        assert fr.total_records == 2 and fr.records_expected == 1
        assert fr.status == ex.FileStatus.PARTIAL
        assert any("mismatch" in r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING)

    def test_intact_file_reconciles_cleanly(self, tmp_path, caplog):
        fp = write(tmp_path / "ok.gb", make_record("TST000026"), make_record("TST000027"))
        out = tmp_path / "out"
        out.mkdir()
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            summary = ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=True)
        fr = summary.file_results[0]
        assert fr.status == ex.FileStatus.SUCCESS and fr.records_expected == 2
        assert not [r for r in caplog.records if r.levelno >= logging.WARNING]

    def test_first_record_broken_is_failed(self, tmp_path):
        fp = write(tmp_path / "bad.gb", make_record("TST000031", length_token="abc"), make_record("TST000032"))
        out = tmp_path / "out"
        out.mkdir()
        summary = ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=True)
        assert summary.file_results[0].status == ex.FileStatus.FAILED

    def test_extract_warns_about_invalid_files_and_always_writes_report(self, tmp_path, caplog):
        d = tmp_path / "in"
        write(d / "ok.gb", make_record("TST000041"))
        (d / "batch_0002.gb").write_text('<?xml version="1.0"?>\n<eSearchResult><ERROR>x</ERROR></eSearchResult>\n')
        out = tmp_path / "out"
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            res = ex.extract([str(d)], str(out), stream=True)
        assert res.success and res.rows == 1
        warned = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
        assert any("batch_0002.gb" in m for m in warned), warned
        report = json.loads((out / "extraction_report.json").read_text(encoding="utf-8"))
        files = {f["file_name"]: f for f in report["files"]}
        assert files["batch_0002.gb"]["status"] == "invalid"
        assert files["batch_0002.gb"]["error_message"]
        assert files["ok.gb"]["status"] == "success"
        assert files["ok.gb"]["total_records"] == 1
        assert files["ok.gb"]["metadata_extracted"] == 1

    def test_extract_reports_partial_file_but_stays_successful(self, tmp_path, caplog):
        d = tmp_path / "in"
        self._broken_file(d)
        out = tmp_path / "out"
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            res = ex.extract([str(d)], str(out), stream=True)
        assert res.success and res.rows == 1
        report = json.loads((out / "extraction_report.json").read_text(encoding="utf-8"))
        (fr,) = report["files"]
        assert fr["status"] == "partial"
        assert fr["record_errors"]
        assert any(r.levelno >= logging.WARNING and "partial" in r.getMessage().lower()
                   for r in caplog.records)

    def test_extract_all_files_failed_returns_failure_but_writes_report(self, tmp_path):
        d = tmp_path / "in"
        (d).mkdir()
        (d / "junk.gb").write_text("not a genbank file\n")
        write(d / "bad.gb", make_record("TST000051", length_token="abc"))
        out = tmp_path / "out"
        res = ex.extract([str(d)], str(out), stream=True)
        assert res.success is False
        report = json.loads((out / "extraction_report.json").read_text(encoding="utf-8"))
        statuses = {f["file_name"]: f["status"] for f in report["files"]}
        assert statuses == {"junk.gb": "invalid", "bad.gb": "failed"}


# ---------------------------------------------------------------------------
# L-04  final.csv must not rewrite numeric-looking text
# ---------------------------------------------------------------------------
class TestL04NumericText:
    def test_combine_preserves_leading_zeros_and_big_numbers(self, tmp_path):
        md = tmp_path / "metadata.csv"
        asm = tmp_path / "assembly.csv"
        out = tmp_path / "final.csv"
        pd.DataFrame({
            "ACCESSION": ["A1", "A2", "A3"],
            "specimen_voucher": ["007", "0123", "45"],
            "strain": ["1E5", "2E5", ""],
            "isolate": ["NA", "", "null"],
            "Ref1PubMed": ["32671913", "", "12345678901234567890"],
        }).to_csv(md, index=False)
        pd.DataFrame({
            "ACCESSION": ["A1", "A3"], "Assembly Name": ["0042", "007"],
        }).to_csv(asm, index=False, sep="\t")
        ex.combine_and_save_final_csv(str(md), str(asm), str(out))
        got = read_str_csv(out)
        assert list(got["specimen_voucher"]) == ["007", "0123", "45"]
        assert list(got["strain"]) == ["1E5", "2E5", ""]
        assert list(got["isolate"]) == ["NA", "", "null"]
        assert list(got["Ref1PubMed"]) == ["32671913", "", "12345678901234567890"]
        assert list(got["Assembly Name"]) == ["0042", "", "007"]

    @pytest.mark.parametrize("stream", [True, False])
    def test_end_to_end_final_csv_keeps_text(self, tmp_path, stream):
        fp = write(
            tmp_path / "num.gb",
            make_record("TST000061", quals={"specimen_voucher": "007", "strain": "1E5"}, pubmed="32671913"),
            make_record("TST000062", quals={"specimen_voucher": "0123", "strain": "2E5"}),
            make_record("TST000063", quals={"specimen_voucher": "45"}),
        )
        out = tmp_path / "out"
        out.mkdir()
        ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=stream)
        got = read_str_csv(out / "final.csv")
        assert list(got["specimen_voucher"]) == ["007", "0123", "45"]
        assert list(got["strain"]) == ["1E5", "2E5", ""]
        assert list(got["Ref1PubMed"]) == ["32671913", "", ""]


# ---------------------------------------------------------------------------
# L-05  taxonomy fallback via Entrez when ete3 is unusable
# ---------------------------------------------------------------------------
PRIAP = {
    "6232": (
        [("phylum", "Priapulida"), ("class", "Priapulimorpha"), ("order", "Priapulimorphida"),
         ("family", "Priapulidae"), ("genus", "Priapulus")],
        ("species", "Priapulus caudatus"),
    ),
}


class TestL05EntrezFallback:
    def _input(self, tmp_path, n_records: int = 2):
        recs = [make_record(f"TST{100 + i:06d}") for i in range(n_records)]
        return write(tmp_path / "in" / "t.gb", *recs)

    def test_fallback_fills_ranks_and_caches(self, tmp_path, monkeypatch):
        fake = install_fake_entrez(monkeypatch, FakeEntrez(PRIAP))
        monkeypatch.setenv("NCBI_EMAIL", "me@example.org")
        monkeypatch.setenv("NCBI_API_KEY", "KEY123")
        fp = self._input(tmp_path)
        out = tmp_path / "out"
        res = ex.extract([fp], str(out), stream=True)
        assert res.success
        got = read_str_csv(out / "final.csv")
        assert set(got["Class"]) == {"Priapulimorpha"}
        assert set(got["Order"]) == {"Priapulimorphida"}
        assert set(got["Family"]) == {"Priapulidae"}
        assert set(got["Genus"]) == {"Priapulus"}
        assert fake.calls == [["6232"]]  # one batched request for the unique TaxonID
        assert fake.seen_email == ["me@example.org"] and fake.seen_api_key == ["KEY123"]
        cache = json.loads((out / "taxonomy_cache.json").read_text(encoding="utf-8"))
        assert cache["6232"]["Genus"] == "Priapulus"

    def test_rerun_uses_cache_without_requests(self, tmp_path, monkeypatch):
        fake = install_fake_entrez(monkeypatch, FakeEntrez(PRIAP))
        fp = self._input(tmp_path)
        out = tmp_path / "out"
        ex.extract([fp], str(out), stream=True)
        assert len(fake.calls) == 1
        # new process state, same output dir -> cache file is reused
        monkeypatch.setattr(ex, "_taxonomy", ex.TaxonomyService())
        res = ex.extract([fp], str(out), stream=True)
        assert res.success
        assert len(fake.calls) == 1, "second run must not hit Entrez again"
        assert set(read_str_csv(out / "final.csv")["Genus"]) == {"Priapulus"}

    def test_requests_are_batched_at_most_200(self, tmp_path, monkeypatch):
        n = 450
        lineages = {
            str(5000 + i): ([("class", "C"), ("genus", f"G{i}")], ("species", f"S{i}"))
            for i in range(n)
        }
        fake = install_fake_entrez(monkeypatch, FakeEntrez(lineages))
        recs = [make_record(f"TST{200 + i:06d}", taxon=str(5000 + i), organism=f"Sp {i}") for i in range(n)]
        fp = write(tmp_path / "in" / "many.gb", *recs)
        out = tmp_path / "out"
        res = ex.extract([fp], str(out), stream=True)
        assert res.success and res.rows == n
        assert [len(c) for c in fake.calls] == [200, 200, 50]
        got = read_str_csv(out / "final.csv")
        assert list(got["Genus"]) == [f"G{i}" for i in range(n)]

    def test_network_failure_warns_and_continues(self, tmp_path, monkeypatch, caplog):
        install_fake_entrez(monkeypatch, FakeEntrez(PRIAP, fail=True))
        fp = self._input(tmp_path)
        out = tmp_path / "out"
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            res = ex.extract([fp], str(out), stream=True)
        assert res.success and res.rows == 2
        got = read_str_csv(out / "final.csv")
        assert set(got["Class"]) == {""} and set(got["Genus"]) == {""}
        assert set(got["TaxonID"]) == {"6232"}
        assert any("Class/Order/Family/Genus will be empty" in r.getMessage() for r in caplog.records)
        # failures are not cached, so the next run retries
        cache_file = out / "taxonomy_cache.json"
        if cache_file.exists():
            assert "6232" not in json.loads(cache_file.read_text(encoding="utf-8"))

    def test_ete3_not_installed_falls_back_without_prompt(self, tmp_path, monkeypatch):
        real_find_spec = importlib.util.find_spec
        monkeypatch.setattr(importlib.util, "find_spec",
                            lambda name, *a, **k: None if name == "ete3" else real_find_spec(name, *a, **k))

        def no_prompt(*a, **k):
            raise AssertionError("must not prompt when the Entrez fallback is available")

        monkeypatch.setattr("builtins.input", no_prompt)
        fake = install_fake_entrez(monkeypatch, FakeEntrez(PRIAP))
        fp = self._input(tmp_path)
        out = tmp_path / "out"
        res = ex.extract([fp], str(out), stream=True)
        assert res.success
        assert fake.calls == [["6232"]]
        assert set(read_str_csv(out / "final.csv")["Family"]) == {"Priapulidae"}


# ---------------------------------------------------------------------------
# L-19  non-taxon db_xref entries are kept
# ---------------------------------------------------------------------------
class TestL19DbXref:
    def test_parse_source_feature_keeps_non_taxon_xrefs(self):
        rec = parse_one(make_record("TST000301", db_xrefs=["taxon:549204", "BOLD:CCANN086-07.COI-5P", "BOLD:XYZ"]))
        src = next(f for f in rec.features if f.type == "source")
        data = ex.parse_source_feature(src)
        assert data["TaxonID"] == "549204"
        assert data["db_xref"] == "BOLD:CCANN086-07.COI-5P; BOLD:XYZ"

    def test_taxon_only_leaves_db_xref_empty(self):
        rec = parse_one(make_record("TST000302", db_xrefs=["taxon:549204"]))
        md = ex.extract_metadata(rec, include_taxonomy=False)
        assert md["TaxonID"] == "549204"
        assert md.get("db_xref", "") == ""

    def test_end_to_end_db_xref_column(self, tmp_path):
        fp = write(
            tmp_path / "x.gb",
            make_record("TST000303", db_xrefs=["taxon:549204", "BOLD:CCANN086-07.COI-5P"]),
            make_record("TST000304"),
        )
        out = tmp_path / "out"
        out.mkdir()
        ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=True)
        got = read_str_csv(out / "final.csv")
        assert list(got["db_xref"]) == ["BOLD:CCANN086-07.COI-5P", ""]
        assert list(got["TaxonID"]) == ["549204", "6232"]


# ---------------------------------------------------------------------------
# L-20  stream and non-stream outputs have identical columns
# ---------------------------------------------------------------------------
class TestL20ColumnParity:
    def test_stream_and_memory_final_csv_have_same_columns_and_values(self, tmp_path):
        fp = write(
            tmp_path / "x.gb",
            make_record("TST000401", quals={"isolate": "I1"}),
            make_record("TST000402", quals={"country": "Russia", "zzz_custom": "q"}),
        )
        frames = {}
        for mode in (True, False):
            out = tmp_path / ("stream" if mode else "memory")
            out.mkdir()
            ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=mode)
            frames[mode] = read_str_csv(out / "final.csv")
        assert list(frames[True].columns) == list(frames[False].columns)
        # columns that no record has are still present (e.g. for downstream --organelle filters)
        for col in ("organelle", "host", "db_xref", "culture_collection", "clone", "Class",
                    "Genus", "Assembly Method", "Expected Final Version"):
            assert col in frames[False].columns
        pd.testing.assert_frame_equal(frames[True], frames[False])

    def test_metadata_csv_columns_match_too(self, tmp_path):
        fp = write(tmp_path / "x.gb", make_record("TST000411", quals={"host": "h"}))
        cols = {}
        for mode in (True, False):
            out = tmp_path / ("s" if mode else "m")
            out.mkdir()
            ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=mode)
            cols[mode] = list(read_str_csv(out / "metadata.csv").columns)
        assert cols[True] == cols[False]


# ---------------------------------------------------------------------------
# L-21  multiple / missing source features
# ---------------------------------------------------------------------------
class TestL21SourceFeature:
    def test_multiple_sources_use_first_and_warn(self, caplog):
        rec = parse_one(make_record(
            "TST000501", quals={"specimen_voucher": "ZMMU:WS30980"},
            extra_sources=[{"organism": "Symbiont sp.", "specimen_voucher": "SYM1", "db_xref": "taxon:999"}]))
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            md = ex.extract_metadata(rec, include_taxonomy=False)
        assert md["specimen_voucher"] == "ZMMU:WS30980"
        assert md["TaxonID"] == "6232"
        assert md["organism"] == "Priapulus caudatus"
        msgs = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
        assert any("TST000501" in m and "source" in m.lower() for m in msgs), msgs

    def test_single_source_does_not_warn(self, caplog):
        rec = parse_one(make_record("TST000502"))
        with caplog.at_level(logging.WARNING, logger=LOGGER_NAME):
            ex.extract_metadata(rec, include_taxonomy=False)
        assert not [r for r in caplog.records if r.levelno >= logging.WARNING]

    def test_missing_source_falls_back_to_annotation_organism(self):
        rec = parse_one(make_record("TST000503", organism="Priapulopsis bicaudatus", source=False))
        assert not [f for f in rec.features if f.type == "source"]
        md = ex.extract_metadata(rec, include_taxonomy=False)
        assert md["Organism"] == "Priapulopsis bicaudatus"
        assert md["organism"] == "Priapulopsis bicaudatus"

    def test_source_organism_wins_over_annotation(self):
        rec = parse_one(make_record("TST000504", quals={"organism": "Source Name"}))
        md = ex.extract_metadata(rec, include_taxonomy=False)
        assert md["organism"] == "Source Name"


# ---------------------------------------------------------------------------
# L-24  assembly_failed is a live counter
# ---------------------------------------------------------------------------
class TestL24AssemblyFailed:
    def test_assembly_failure_is_counted_and_summarised(self, tmp_path, caplog):
        fp = write(tmp_path / "a.gb", make_record("TST000601", assembly={"Assembly Name": "asm1"}))

        def boom(md, asm, src, has, track, result):
            raise RuntimeError("disk full while writing assembly row")

        fr = ex._process_record_loop(fp, boom, False, True, True)
        assert fr.assembly_failed == 1
        assert fr.assembly_extracted == 0
        assert fr.metadata_failed == 1

        summary = ex.ProcessingSummary()
        summary.add_result(fr)
        assert summary.total_assembly_failed == 1
        with caplog.at_level(logging.INFO, logger=LOGGER_NAME):
            summary.print_summary()
        text = "\n".join(r.getMessage() for r in caplog.records)
        assert "Assembly coverage:     0/1 (0.0%)" in text

    def test_successful_assembly_is_not_counted_as_failed(self, tmp_path):
        fp = write(tmp_path / "a.gb", make_record("TST000611", assembly={"Assembly Name": "asm1"}))
        out = tmp_path / "out"
        out.mkdir()
        summary = ex.process_all_files([fp], str(out), include_taxonomy=False, stream_mode=False)
        assert summary.total_assembly_extracted == 1
        assert summary.total_assembly_failed == 0
