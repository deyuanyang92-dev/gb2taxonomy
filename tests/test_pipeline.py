"""Tests for _pipeline.py: PipelineResult and _find_fallback."""

import pytest
from pathlib import Path

from g2t._pipeline import PipelineResult, _find_fallback


class TestPipelineResult:
    def test_creation(self):
        r = PipelineResult(success=True, output_dir="/tmp/out")
        assert r.success is True
        assert r.output_dir == "/tmp/out"
        assert r.steps_completed == 0
        assert r.final_output == ""
        assert r.elapsed == 0.0
        assert r.log == []

    def test_with_all_fields(self):
        r = PipelineResult(
            success=False,
            output_dir="/tmp/out",
            steps_completed=2,
            final_output="/tmp/out/organized.csv",
            elapsed=5.3,
            log=["Step 1 done", "Step 2 failed"],
        )
        assert r.success is False
        assert r.steps_completed == 2
        assert len(r.log) == 2

    def test_log_mutable_default(self):
        r1 = PipelineResult(success=True, output_dir="/a")
        r2 = PipelineResult(success=True, output_dir="/b")
        r1.log.append("msg")
        assert r2.log == []  # Separate instances have separate lists


class TestFindFallback:
    def test_raises_when_no_candidates(self, tmp_path):
        target = tmp_path / "out.csv"
        with pytest.raises(FileNotFoundError, match="no valid output"):
            _find_fallback(tmp_path, "*.csv", target, "TestStep")

    def test_no_error_when_target_valid(self, tmp_path):
        target = tmp_path / "out.csv"
        target.write_text("col\nval\n")
        _find_fallback(tmp_path, "*.csv", target, "TestStep")

    def test_raises_when_only_empty_files(self, tmp_path):
        empty = tmp_path / "empty.csv"
        empty.write_text("")
        target = tmp_path / "out.csv"
        with pytest.raises(FileNotFoundError):
            _find_fallback(tmp_path, "*.csv", target, "TestStep")
