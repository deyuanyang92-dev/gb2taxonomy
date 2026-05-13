"""Tests for custom exceptions and sys.exit removal (Phase 3)."""

import pytest

from g2t.utils import G2TError, UnsupportedFormatError, PipelineStepError


class TestExceptionHierarchy:
    def test_unsupported_format_is_g2t_error(self):
        with pytest.raises(G2TError):
            raise UnsupportedFormatError("test")

    def test_pipeline_step_is_g2t_error(self):
        with pytest.raises(G2TError):
            raise PipelineStepError("test")

    def test_unsupported_format_message(self):
        exc = UnsupportedFormatError(".json")
        assert ".json" in str(exc)


class TestClassifyNoSysExit:
    def test_process_prematch_unsupported_format_raises(self, tmp_path):
        """Verify classify.py raises UnsupportedFormatError instead of sys.exit."""
        from g2t.classify import process_prematch, MatchConfig

        json_file = tmp_path / "test.json"
        json_file.write_text('{"a": 1}')

        with pytest.raises(UnsupportedFormatError):
            process_prematch(
                input_file=str(json_file),
                output_dir=str(tmp_path / "out"),
                config=MatchConfig(),
                prematch_output="prematch.csv",
                multip_output="multip.csv",
                unmatched_output="unmatched.csv",
                conflicted_output="conflicted.csv",
                moleculetype="all",
                mol_type="all",
                length_filter="all",
                organelle_filter="all",
            )

    def test_process_recheck_unsupported_format_raises(self, tmp_path):
        """Verify recheck raises UnsupportedFormatError instead of sys.exit."""
        from g2t.classify import process_recheck, MatchConfig

        json_file = tmp_path / "test.json"
        json_file.write_text('{"a": 1}')

        with pytest.raises(UnsupportedFormatError):
            process_recheck(
                input_file=str(json_file),
                output_dir=str(tmp_path / "out"),
                config=MatchConfig(),
            )
