"""
Accession-level store of downloaded GenBank records
====================================================
One SQLite file maps ``accession.version`` -> flat-file record text (zlib-compressed).
``g2t.download`` fetches only the accession.versions that are not in the store yet, so
re-running a download, choosing another selection of the same taxon, or downloading a
taxon that overlaps an earlier one (e.g. a family inside an order) reuses every record
already fetched. A new version of a record (``AB123456.2``) is a new key and is fetched.
"""

from __future__ import annotations

import re
import sqlite3
import zlib
from collections.abc import Iterable, Iterator
from datetime import date
from pathlib import Path

_VERSION = re.compile(r"^VERSION\s+(\S+)")


def iter_records(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    """Split GenBank flat-file lines into (accession.version, record text); records end with '//'."""
    buf: list[str] = []
    acc = ""
    for line in lines:
        if line.startswith("LOCUS"):
            buf, acc = [], ""
        buf.append(line)
        if not acc:
            m = _VERSION.match(line)
            if m:
                acc = m.group(1)
        if line.startswith("//"):
            if acc:
                text = "".join(buf)
                yield acc, text if text.endswith("\n") else text + "\n"
            buf, acc = [], ""


class RecordStore:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(str(self.path))
        self.db.execute("CREATE TABLE IF NOT EXISTS records (acc TEXT PRIMARY KEY, fetched TEXT, data BLOB)")
        self.db.commit()

    def __contains__(self, acc: str) -> bool:
        return self.db.execute("SELECT 1 FROM records WHERE acc=?", (acc,)).fetchone() is not None

    def missing(self, accs: Iterable[str]) -> list[str]:
        have: set[str] = set()
        lst = list(accs)
        for i in range(0, len(lst), 500):
            chunk = lst[i:i + 500]
            q = "SELECT acc FROM records WHERE acc IN ({})".format(",".join("?" * len(chunk)))
            have.update(r[0] for r in self.db.execute(q, chunk))
        return [a for a in lst if a not in have]

    def get(self, acc: str) -> str | None:
        r = self.db.execute("SELECT data FROM records WHERE acc=?", (acc,)).fetchone()
        return zlib.decompress(r[0]).decode("utf-8") if r else None

    def put_many(self, items: Iterable[tuple[str, str]]) -> int:
        today = date.today().isoformat()
        rows = [(a, today, zlib.compress(t.encode("utf-8"), 6)) for a, t in items]
        self.db.executemany("INSERT OR REPLACE INTO records VALUES (?,?,?)", rows)
        self.db.commit()
        return len(rows)

    def import_files(self, files: Iterable[Path], wanted: set[str]) -> int:
        """Copy records listed in ``wanted`` from existing batch files (streamed, so huge files are fine)."""
        n = 0
        for f in files:
            if not wanted:
                break
            pending: list[tuple[str, str]] = []
            with open(f, encoding="utf-8", errors="ignore") as fh:
                for acc, text in iter_records(fh):
                    if acc in wanted:
                        pending.append((acc, text))
                        wanted.discard(acc)
                    if len(pending) >= 500:
                        n += self.put_many(pending)
                        pending = []
            n += self.put_many(pending)
        return n

    def close(self) -> None:
        self.db.close()
