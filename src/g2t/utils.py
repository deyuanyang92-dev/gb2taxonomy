"""Shared utilities for the g2t pipeline."""

from __future__ import annotations

import argparse
import os
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd


@dataclass
class StepResult:
    """Standardized result object returned by each pipeline step."""
    success: bool
    output_file: str
    rows: int = 0
    elapsed: float = 0.0


class G2TError(Exception):
    """Base exception for g2t pipeline."""
    pass


class UnsupportedFormatError(G2TError):
    """Raised when input file format is not supported."""
    pass


class PipelineStepError(G2TError):
    """Raised when a pipeline step encounters an unrecoverable error."""
    pass


def str2bool(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in ("yes", "true", "t", "1", "y", "on"):
        return True
    if s in ("no", "false", "f", "0", "n", "off"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got: {x!r}")


def sanitize_string(s: Any) -> str:
    if s is None:
        return ""
    try:
        if pd.isna(s):
            return ""
    except (ValueError, TypeError):
        pass
    s = str(s).strip()
    if not s:
        return ""
    s = unicodedata.normalize("NFKC", s)
    s = s.replace("'", "").replace('"', "")
    s = s.replace("(", "").replace(")", "")
    s = s.replace("　", " ")
    s = s.replace(" ", "_")
    s = re.sub(r"[<>:./\\|?*]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def safe_concat(a: Any, b: Any) -> str:
    sa = sanitize_string(a)
    sb = sanitize_string(b)
    if sa and sb:
        return f"{sa}_{sb}"
    return sa or sb


def read_table(file_path: str) -> pd.DataFrame:
    ext = os.path.splitext(file_path)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(file_path, dtype=str)
    elif ext in (".tsv", ".txt"):
        return pd.read_csv(file_path, sep="\t", dtype=str)
    else:
        return pd.read_csv(file_path, dtype=str)


def write_csv(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8")


def parse_interval(interval_str: str) -> Tuple[Optional[int], Optional[int]]:
    if not interval_str:
        return None, None
    sep = ":" if ":" in interval_str else ","
    parts = interval_str.split(sep, 1)
    lo = int(parts[0]) if parts[0].strip() else None
    hi = int(parts[1]) if len(parts) > 1 and parts[1].strip() else None
    return lo, hi


def length_in_range(length_val: Any, range_str: str) -> bool:
    if not range_str or range_str.lower() == "none":
        return True
    lo, hi = parse_interval(range_str)
    try:
        val = int(length_val)
    except (ValueError, TypeError):
        return True
    if lo is not None and val < lo:
        return False
    if hi is not None and val > hi:
        return False
    return True


def is_valid_output(path: Union[str, Path]) -> bool:
    p = Path(path)
    return p.exists() and p.is_file() and p.stat().st_size > 0


def col_key(s: str) -> str:
    return s.lower().replace(" ", "_")


def resolve_col(df: pd.DataFrame, name: str, lookup: Optional[Dict[str, str]] = None) -> Optional[str]:
    if name in df.columns:
        return name
    if lookup is None:
        lookup = {col_key(c): c for c in df.columns}
    k = col_key(name)
    if k in lookup:
        return lookup[k]
    alt = name.replace("_", " ")
    if alt in df.columns:
        return alt
    k2 = col_key(alt)
    return lookup.get(k2)


def build_col_lookup(df: pd.DataFrame) -> Dict[str, str]:
    return {col_key(c): c for c in df.columns}
