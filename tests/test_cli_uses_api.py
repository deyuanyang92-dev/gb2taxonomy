"""g2t-extract and g2t-classify run the same functions as the pipeline (one code path each)."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest
from Bio import Entrez

from g2t.utils import StepResult

ex = importlib.import_module("g2t.extract")
_tax = importlib.import_module("g2t._extract.taxonomy")  # the taxonomy service lives here
cl = importlib.import_module("g2t.classify")

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "haplotaxidae" / "haplotaxidae.gb"


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    monkeypatch.setattr(_tax, "_taxonomy", ex.TaxonomyService())
    monkeypatch.setitem(sys.modules, "ete3", None)

    def blocked(*_a, **_k):
        raise OSError("network access is blocked in tests")

    monkeypatch.setattr(Entrez, "efetch", blocked)
    monkeypatch.setattr(Entrez, "sleep_between_tries", 0)


@pytest.fixture
def one_record(tmp_path):
    text = EXAMPLE.read_text(encoding="utf-8")
    f = tmp_path / "one.gb"
    f.write_text(text[: text.index("\n//") + 4], encoding="utf-8")
    return f


def _capture(monkeypatch, module, name):
    seen = {}

    def fake(*args, **kwargs):
        seen["args"], seen["kwargs"] = args, kwargs
        return StepResult(success=True, output_file="x", rows=1)

    monkeypatch.setattr(module, name, fake)
    return seen


def test_extract_cli_calls_extract(monkeypatch, tmp_path):
    seen = _capture(monkeypatch, ex, "extract")
    ex.main(["-i", "a.gb", "-o", str(tmp_path), "-oam", "--no_tracking", "--final_delimiter", ";"])
    assert seen["args"] == (["a.gb"], str(tmp_path))
    kw = seen["kwargs"]
    assert kw["only_assembly"] and not kw["only_metadata"] and kw["details"]
    assert not kw["tracking"] and kw["report"] and kw["final_delimiter"] == ";" and kw["max_tasks"] == 4


def test_classify_cli_calls_classify(monkeypatch, tmp_path):
    seen = _capture(monkeypatch, cl, "classify")
    cl.main(["-i", "in.csv", "-o", str(tmp_path), "--organelle", "mitochondrion", "-wgte", "COI,16s"])
    kw = seen["kwargs"]
    assert seen["args"][:2] == ("in.csv", str(tmp_path)) and isinstance(seen["args"][2], cl.MatchConfig)
    assert kw["organelle"] == "mitochondrion" and kw["extract_gene_types"] == ["coi", "16s"]


@pytest.mark.parametrize("module", [ex, cl])
def test_cli_exits_nonzero_when_the_step_fails(monkeypatch, tmp_path, module):
    name = "extract" if module is ex else "classify"
    monkeypatch.setattr(module, name, lambda *_a, **_k: StepResult(success=False, output_file=""))
    with pytest.raises(SystemExit) as e:
        module.main(["-i", "x", "-o", str(tmp_path)])
    assert e.value.code == 1


def test_extract_cli_fails_when_no_genbank_file(tmp_path):
    (tmp_path / "notes.txt").write_text("hello\n", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        ex.main(["-i", str(tmp_path), "-o", str(tmp_path / "out")])
    assert e.value.code == 1


def test_classify_failure_is_logged_with_traceback(tmp_path, caplog):
    r = cl.classify(str(tmp_path / "missing.csv"), str(tmp_path / "out"))
    assert not r.success
    assert any(rec.exc_info for rec in caplog.records if "Classify failed" in rec.getMessage())


@pytest.mark.parametrize(
    "mode,product,rows",
    # this record has no assembly comment, so the assembly-only table has no row for it
    [({"only_assembly": True}, "assembly.csv", 0), ({"only_metadata": True}, "metadata.csv", 1)],
)
def test_single_half_modes_succeed_with_their_own_table(tmp_path, one_record, mode, product, rows):
    """-oam / -oeam write no merged table; their own table is the product."""
    r = ex.extract([str(one_record)], str(tmp_path), **mode)
    assert r.success and r.output_file.endswith(product) and r.rows == rows
    assert not (tmp_path / "final.csv").exists()


def test_cli_writes_its_extra_reports(tmp_path, one_record):
    ex.main(["-i", str(one_record), "-o", str(tmp_path)])
    assert (tmp_path / "final.csv").exists() and (tmp_path / "record_tracking.csv").exists()
    assert (tmp_path / "extraction_report.json").exists()


def test_pipeline_call_writes_no_cli_extras(tmp_path, one_record):
    r = ex.extract([str(one_record)], str(tmp_path))
    assert r.success and r.rows == 1
    assert not (tmp_path / "record_tracking.csv").exists()


def test_python_m_classify_writes_classify_log(tmp_path):
    """Under ``python -m g2t.classify`` the module is ``__main__``; its log lines still reach classify.log."""
    import subprocess

    import pandas as pd

    f = tmp_path / "final.csv"
    pd.DataFrame(
        [
            {
                "LocusID": "A",
                "ACCESSION": "A",
                "Length": 658,
                "Definition": "Priapulus caudatus cytochrome c oxidase subunit I (COI) gene; mitochondrial",
                "Organism": "Priapulus caudatus",
            }
        ]
    ).to_csv(f, index=False)
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, "-m", "g2t.classify", "-i", str(f), "-o", str(out)], check=True, capture_output=True
    )
    log = (out / "classify.log").read_text(encoding="utf-8")
    assert "Gene Type Classifier" in log and "Round 1: Gene type matching" in log


def test_pipeline_failure_says_why(tmp_path):
    """The pipeline result names the cause, not only 'Step 1 failed'."""
    from g2t._pipeline import run

    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "notes.txt").write_text("hello\n", encoding="utf-8")
    r = run([str(tmp_path / "in")], str(tmp_path / "out"), quiet=True)
    assert not r.success
    assert "Step 1 failed: no GenBank files" in r.log[-1]


def test_classify_failure_carries_the_error(tmp_path):
    r = cl.classify(str(tmp_path / "missing.csv"), str(tmp_path / "out"))
    assert not r.success and r.error


def test_organize_failure_names_missing_columns(tmp_path):
    import pandas as pd

    from g2t.organize import organize

    f = tmp_path / "in.csv"
    pd.DataFrame({"LocusID": ["A"]}).to_csv(f, index=False)
    r = organize(str(f), str(tmp_path / "out.csv"))
    assert not r.success and "missing required columns" in r.error and "gene_type" in r.error
