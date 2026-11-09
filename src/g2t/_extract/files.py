"""Finding, opening and checking GenBank input files."""

from __future__ import annotations

import logging
import os
import shutil

logger = logging.getLogger("g2t.extract")  # one logger name for the whole step


try:
    import psutil
except ImportError:
    psutil = None

GB_EXTENSIONS = (".gb", ".gbf", ".gbk", ".gbff", ".genbank")


def is_genbank_name(name: str) -> bool:
    """A GenBank file name: one of GB_EXTENSIONS, optionally gzip-compressed (``.gb.gz``), case-insensitive."""
    n = name.lower()
    return n.endswith(GB_EXTENSIONS) or n.endswith(tuple(e + ".gz" for e in GB_EXTENSIONS))


def open_genbank(path: str, binary: bool = False):
    """Open a GenBank file for reading, decompressing ``.gz`` files on the fly."""
    if path.lower().endswith(".gz"):
        import gzip

        return gzip.open(path, "rb") if binary else gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rb") if binary else open(path, encoding="utf-8", errors="replace")  # noqa: SIM115


# ===============================
# File discovery and validation
# ===============================
def discover_gb_files(input_paths: list[str], recursive: bool = True) -> tuple[list[str], int, int]:
    gb_files = []
    total_files = 0
    skipped = 0
    for path in input_paths:
        if os.path.isdir(path):
            walker = (
                os.walk(path)
                if recursive
                else [(path, [], [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])]
            )
            for root, dirs, files in walker:
                dirs.sort()  # os.walk honours in-place sorting: same file order on every filesystem
                for fn in sorted(files):
                    total_files += 1
                    if is_genbank_name(fn):
                        gb_files.append(os.path.join(root, fn))
                    else:
                        skipped += 1
        elif os.path.isfile(path):
            total_files += 1
            if is_genbank_name(path):
                gb_files.append(path)
            else:
                skipped += 1
        else:
            logger.warning(f"Path not found: {path}")
    logger.info(f"File discovery: {total_files} files, {len(gb_files)} GenBank, {skipped} skipped")
    return gb_files, total_files, skipped


# Combined validation + record counting in single file read
def validate_and_estimate_records(file_path: str) -> tuple[bool, str, int]:
    """Single-pass: validate file header AND count '//' delimiters."""
    try:
        if not os.path.exists(file_path):
            return False, "File does not exist", 0
        sz = os.path.getsize(file_path)
        if sz == 0:
            return False, "Empty file (0 bytes)", 0

        count = 0
        with open_genbank(file_path) as f:
            header = f.read(2048)
            if not header.strip():
                return False, "File content empty", 0
            if not header.lstrip().startswith("LOCUS"):
                return False, "Does not start with LOCUS", 0
            # Continue reading for record count
            f.seek(0)
            for line in f:
                if line.strip() == "//":
                    count += 1
        return True, "", count
    except PermissionError:
        return False, "No read permission", 0
    except Exception as e:
        return False, f"Validation error: {e}", 0


def format_size(size_bytes) -> str:
    s = float(size_bytes) if size_bytes else 0
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if s < 1024:
            return f"{s:.1f} {unit}"
        s /= 1024
    return f"{s:.1f} PB"


def copy_failed_files(failed_files: list[str], output_dir: str):
    if not failed_files:
        return
    failed_dir = os.path.join(output_dir, "failed_extract")
    os.makedirs(failed_dir, exist_ok=True)
    copied = 0
    for fp in failed_files:
        try:
            dest = os.path.join(failed_dir, os.path.basename(fp))
            if os.path.exists(dest):
                base, ext = os.path.splitext(os.path.basename(fp))
                c = 1
                while os.path.exists(dest):
                    dest = os.path.join(failed_dir, f"{base}_{c}{ext}")
                    c += 1
            shutil.copy2(fp, dest)
            copied += 1
        except OSError as e:
            logger.error(f"Failed to copy {fp}: {e}")
    logger.info(f"Archived {copied}/{len(failed_files)} failed files -> {failed_dir}")


def get_memory_mb() -> float | None:
    if psutil:
        try:
            return float(psutil.Process(os.getpid()).memory_info().rss) / (1024 * 1024)
        except Exception:
            pass
    return None
