"""End-to-end run on the bundled example (examples/haplotaxidae, 347 GenBank records), offline."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest
from Bio import Entrez

import g2t
from g2t import _cli
from g2t import extract as ex
from g2t._extract import taxonomy as _tax
from g2t._pipeline import run

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "haplotaxidae" / "haplotaxidae.gb"


def _blocked(*args, **kwargs):
    raise OSError("network access is blocked in tests")


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    mp.setattr(_tax, "_taxonomy", ex.TaxonomyService())
    mp.setitem(sys.modules, "ete3", None)
    mp.setattr(Entrez, "efetch", _blocked)
    mp.setattr(Entrez, "sleep_between_tries", 0)
    out = tmp_path_factory.mktemp("example")
    res = run(input_files=[str(EXAMPLE)], output_dir=str(out), quiet=True)
    mp.undo()
    return res, out


def test_pipeline_succeeds(result):
    res, out = result
    assert res.success, res.log
    assert res.steps_completed == 6  # 1, 2, 3, 3b, 4, 5 (2b is a check of step 2, not counted)
    assert (out / "organized_genes" / "screen.csv").exists()
    m = pd.read_csv(out / "organized_genes" / "organized_species_voucher.csv")
    assert "screen_flags" in m.columns


def test_every_record_has_a_status(result):
    _, out = result
    st = pd.read_csv(out / "labeled_genes" / "record_status.csv")
    assert len(st) == 347
    assert (st.status == "assigned").all()


def test_specimen_rows_link_genes(result):
    _, out = result
    m = pd.read_csv(out / "organized_genes" / "organized_species_voucher.csv", dtype=str)
    # Haplotaxis sp. SE, specimen CE1233 (Martin et al. 2024, Table S1): five genes in one row
    row = m[m["16s"].fillna("").str.contains("PP994580")]
    assert len(row) == 1
    for col, acc in [("12s", "PP968574"), ("18s", "PP990634"), ("28s", "PP994523"), ("coi", "PP988405")]:
        assert acc in row.iloc[0][col]


def test_run_record_written(result):
    _, out = result
    prov = json.loads((out / "g2t_run.json").read_text())
    assert prov["versions"]["g2t"] == g2t.__version__
    assert prov["parameters"]["input_files"] == [str(EXAMPLE)]
    assert prov["inputs"][0]["bytes"] == EXAMPLE.stat().st_size
    import hashlib

    assert prov["inputs"][0]["sha256"] == hashlib.sha256(EXAMPLE.read_bytes()).hexdigest()
    assert "taxonomy_source" in prov
    assert prov["steps_completed"] == 6


@pytest.mark.parametrize(
    "entry,prog",
    [
        ("main_pipeline", "g2t"),
        ("main_download", "g2t-download"),
        ("main_extract", "g2t-extract"),
        ("main_classify", "g2t-classify"),
        ("main_voucher", "g2t-voucher"),
        ("main_reconcile", "g2t-reconcile"),
        ("main_curate", "g2t-curate"),
        ("main_organize", "g2t-organize"),
    ],
)
def test_version_flag(entry, prog, capsys):
    with pytest.raises(SystemExit) as e:
        getattr(_cli, entry)(["--version"])
    assert e.value.code == 0
    assert capsys.readouterr().out.strip() == f"{prog} {g2t.RELEASE}" == f"{prog} v0.05"


def test_gzipped_input_gives_the_same_matrix(result, tmp_path):
    """A .gb.gz copy of the example gives the same specimen matrix as the plain file."""
    import gzip
    import shutil

    gz = tmp_path / "haplotaxidae.gb.gz"
    with open(EXAMPLE, "rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    mp = pytest.MonkeyPatch()
    mp.setattr(_tax, "_taxonomy", ex.TaxonomyService())
    mp.setitem(sys.modules, "ete3", None)
    mp.setattr(Entrez, "efetch", _blocked)
    mp.setattr(Entrez, "sleep_between_tries", 0)
    res = run(input_files=[str(gz)], output_dir=str(tmp_path / "out"), quiet=True)
    mp.undo()
    assert res.success, res.log
    name = "organized_genes/organized_species_voucher.csv"
    plain = pd.read_csv(result[1] / name, dtype=str)
    zipped = pd.read_csv(tmp_path / "out" / name, dtype=str)
    pd.testing.assert_frame_equal(plain, zipped)
