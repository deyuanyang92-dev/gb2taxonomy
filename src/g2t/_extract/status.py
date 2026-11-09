"""Per-record and per-file status, the run summary and its reports."""

from __future__ import annotations

import csv
import json
import logging
import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

logger = logging.getLogger("g2t.extract")  # one logger name for the whole step


# ===============================
# Data classes
# ===============================
class FileStatus(Enum):
    SUCCESS = "success"
    PARTIAL = "partial"
    FAILED = "failed"
    SKIPPED = "skipped"
    INVALID = "invalid"
    EMPTY = "empty"


@dataclass
class RecordError:
    file_path: str
    record_index: int
    accession: str = ""
    error_type: str = ""
    error_message: str = ""
    phase: str = ""


@dataclass
class RecordTrack:
    accession: str
    record_index: int
    metadata_ok: bool = False
    assembly_ok: bool = False
    assembly_source: str = ""
    has_assembly_data: bool = False
    error_message: str = ""


@dataclass
class FileResult:
    file_path: str
    status: FileStatus = FileStatus.FAILED
    total_records: int = 0
    metadata_extracted: int = 0
    metadata_failed: int = 0
    assembly_extracted: int = 0
    assembly_no_data: int = 0
    assembly_failed: int = 0
    # Number of "//" record terminators found in the file (None = not counted).
    records_expected: int | None = None
    error_message: str = ""
    elapsed_time: float = 0.0
    record_errors: list[RecordError] = field(default_factory=list)
    record_tracks: list[RecordTrack] = field(default_factory=list)

    def add_record_error(self, error: RecordError):
        self.record_errors.append(error)

    @property
    def extracted_records(self) -> int:
        return self.metadata_extracted

    @property
    def failed_records(self) -> int:
        return self.metadata_failed


@dataclass
class ProcessingSummary:
    total_files_found: int = 0
    gb_files_found: int = 0
    non_gb_files_skipped: int = 0
    files_succeeded: int = 0
    files_partial: int = 0
    files_failed: int = 0
    files_empty: int = 0
    files_invalid: int = 0
    total_records_in_files: int = 0
    total_metadata_extracted: int = 0
    total_metadata_failed: int = 0
    total_assembly_extracted: int = 0
    total_assembly_no_data: int = 0
    total_assembly_failed: int = 0
    total_elapsed_time: float = 0.0
    file_results: list[FileResult] = field(default_factory=list)

    def add_result(self, result: FileResult):
        self.file_results.append(result)
        _map = {
            FileStatus.SUCCESS: "files_succeeded",
            FileStatus.PARTIAL: "files_partial",
            FileStatus.FAILED: "files_failed",
            FileStatus.EMPTY: "files_empty",
            FileStatus.INVALID: "files_invalid",
        }
        attr = _map.get(result.status)
        if attr:
            setattr(self, attr, getattr(self, attr) + 1)
        self.total_records_in_files += result.total_records
        self.total_metadata_extracted += result.metadata_extracted
        self.total_metadata_failed += result.metadata_failed
        self.total_assembly_extracted += result.assembly_extracted
        self.total_assembly_no_data += result.assembly_no_data
        self.total_assembly_failed += result.assembly_failed

    def get_failed_files(self) -> list[str]:
        return [
            r.file_path
            for r in self.file_results
            if r.status in (FileStatus.FAILED, FileStatus.INVALID, FileStatus.EMPTY)
        ]

    def get_all_record_errors(self) -> list[RecordError]:
        errs = []
        for r in self.file_results:
            errs.extend(r.record_errors)
        return errs

    def print_summary(self):
        logger.info("")
        logger.info("=" * 75)
        logger.info("                          Summary")
        logger.info("=" * 75)
        logger.info(f"  Total files:              {self.total_files_found}")
        logger.info(f"    GenBank files:          {self.gb_files_found}")
        logger.info(f"    Non-GB (skipped):       {self.non_gb_files_skipped}")
        logger.info("-" * 75)
        logger.info("  File results:")
        logger.info(f"    Success:               {self.files_succeeded}")
        logger.info(f"    Partial:               {self.files_partial}")
        logger.info(f"    Failed:                {self.files_failed}")
        logger.info(f"    Empty:                 {self.files_empty}")
        logger.info(f"    Invalid:               {self.files_invalid}")
        logger.info("-" * 75)
        logger.info("  Record statistics:")
        logger.info(f"    Total records:         {self.total_records_in_files}")
        logger.info(f"    Metadata extracted:    {self.total_metadata_extracted}")
        logger.info(f"    Metadata failed:       {self.total_metadata_failed}")
        md_rate = (
            self.total_metadata_extracted / self.total_records_in_files * 100 if self.total_records_in_files > 0 else 0
        )
        logger.info(f"    Metadata rate:         {md_rate:.1f}%")
        logger.info("")
        logger.info(f"    Assembly extracted:    {self.total_assembly_extracted}")
        logger.info(f"    Assembly no data:      {self.total_assembly_no_data}")
        logger.info(f"    Assembly failed:       {self.total_assembly_failed}")
        asm_with_data = self.total_assembly_extracted + self.total_assembly_failed
        asm_rate = self.total_assembly_extracted / asm_with_data * 100 if asm_with_data > 0 else 0
        logger.info(f"    Assembly coverage:     {self.total_assembly_extracted}/{asm_with_data} ({asm_rate:.1f}%)")
        logger.info("-" * 75)
        logger.info(f"  Total time:              {self.total_elapsed_time:.2f}s")
        logger.info("=" * 75)
        self._print_per_file_details()

    def _print_per_file_details(self):
        if not self.file_results:
            return
        logger.info("")
        logger.info("  Per-file details:")
        logger.info("  " + "-" * 73)
        for r in self.file_results:
            fname = os.path.basename(r.file_path)
            icon = {"success": "OK", "partial": "~", "failed": "X", "empty": "0", "invalid": "!"}.get(
                r.status.value, "?"
            )
            logger.info(f"  {icon} {fname}")
            if r.status in (FileStatus.INVALID, FileStatus.EMPTY):
                logger.info(f"      Status: {r.status.value}  |  Reason: {r.error_message}")
                continue
            md_str = f"Metadata: {r.metadata_extracted}/{r.total_records}"
            asm_str = f"Assembly: {r.assembly_extracted}"
            if r.assembly_no_data > 0:
                asm_str += f" (no_data:{r.assembly_no_data})"
            if r.assembly_failed > 0:
                asm_str += f" (failed:{r.assembly_failed})"
            logger.info(f"      Records: {r.total_records}  |  {md_str}  |  {asm_str}  |  {r.elapsed_time:.1f}s")
            if r.status != FileStatus.SUCCESS and r.error_message:
                logger.info(f"      Note: {r.error_message}")

            if r.metadata_failed > 0 or r.assembly_failed > 0:
                failed_accs = [t.accession for t in r.record_tracks if not t.metadata_ok]
                if failed_accs:
                    shown = failed_accs[:5]
                    more = f" ... +{len(failed_accs) - 5}" if len(failed_accs) > 5 else ""
                    logger.info(f"      Failed ACCESSIONs: {', '.join(shown)}{more}")

    def save_report(self, output_dir: str):
        report_path = os.path.join(output_dir, "extraction_report.json")
        report: dict[str, Any] = {
            "summary": {
                "total_files_found": self.total_files_found,
                "gb_files_found": self.gb_files_found,
                "non_gb_files_skipped": self.non_gb_files_skipped,
                "files_succeeded": self.files_succeeded,
                "files_partial": self.files_partial,
                "files_failed": self.files_failed,
                "files_empty": self.files_empty,
                "files_invalid": self.files_invalid,
                "total_records_in_files": self.total_records_in_files,
                "total_metadata_extracted": self.total_metadata_extracted,
                "total_metadata_failed": self.total_metadata_failed,
                "total_assembly_extracted": self.total_assembly_extracted,
                "total_assembly_no_data": self.total_assembly_no_data,
                "total_assembly_failed": self.total_assembly_failed,
                "total_elapsed_time": round(self.total_elapsed_time, 2),
            },
            "files": [],
        }
        for r in self.file_results:
            fr = {
                "file_path": r.file_path,
                "file_name": os.path.basename(r.file_path),
                "status": r.status.value,
                "total_records": r.total_records,
                "metadata_extracted": r.metadata_extracted,
                "metadata_failed": r.metadata_failed,
                "assembly_extracted": r.assembly_extracted,
                "assembly_no_data": r.assembly_no_data,
                "assembly_failed": r.assembly_failed,
                "records_expected": r.records_expected,
                "error_message": r.error_message,
                "elapsed_time": round(r.elapsed_time, 2),
                "record_errors": [
                    {
                        "record_index": e.record_index,
                        "accession": e.accession,
                        "error_type": e.error_type,
                        "error_message": e.error_message,
                        "phase": e.phase,
                    }
                    for e in r.record_errors
                ],
            }
            report["files"].append(fr)
        try:
            with open(report_path, "w", encoding="utf-8") as f:
                json.dump(report, f, ensure_ascii=False, indent=2)
            logger.info(f"Report: {report_path}")
        except Exception as e:
            logger.error(f"Failed to save report: {e}")

    def save_error_csv(self, output_dir: str):
        all_errors = self.get_all_record_errors()
        if not all_errors:
            return
        path = os.path.join(output_dir, "record_errors.csv")
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["file_name", "record_index", "accession", "error_type", "error_message", "phase"])
                for e in all_errors:
                    w.writerow(
                        [
                            os.path.basename(e.file_path),
                            e.record_index,
                            e.accession,
                            e.error_type,
                            e.error_message,
                            e.phase,
                        ]
                    )
            logger.info(f"Error CSV: {path} ({len(all_errors)} records)")
        except Exception as e:
            logger.error(f"Failed to save error CSV: {e}")

    def save_tracking_csv(self, output_dir: str):
        all_tracks = []
        for r in self.file_results:
            all_tracks.extend(r.record_tracks)
        if not all_tracks:
            return
        path = os.path.join(output_dir, "record_tracking.csv")
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(
                    [
                        "file_name",
                        "record_index",
                        "accession",
                        "metadata_ok",
                        "assembly_ok",
                        "assembly_source",
                        "has_assembly_data",
                        "error_message",
                    ]
                )
                for r in self.file_results:
                    fname = os.path.basename(r.file_path)
                    for t in r.record_tracks:
                        w.writerow(
                            [
                                fname,
                                t.record_index,
                                t.accession,
                                t.metadata_ok,
                                t.assembly_ok,
                                t.assembly_source,
                                t.has_assembly_data,
                                t.error_message,
                            ]
                        )
            logger.info(f"Tracking CSV: {path} ({len(all_tracks)} records)")
        except Exception as e:
            logger.error(f"Failed to save tracking CSV: {e}")
