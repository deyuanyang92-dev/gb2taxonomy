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


class DirStore:
    """Store-free mode: the batch files of one selection directory are the only copy of the records.

    Present records are found by scanning VERSION lines; fetched records are written as new batch files;
    superseded versions and records no longer in the selection are removed from the old files (``prune``).
    """

    def __init__(self, out: str | Path):
        self.out = Path(out)
        self.where: dict[str, Path] = {}
        for f in sorted(self.out.glob("batch_*.gb")):
            with open(f, encoding="utf-8", errors="ignore") as fh:
                for acc, _ in iter_records(fh):
                    self.where[acc] = f
        nums = [int(m.group(1)) for f in self.out.glob("batch_*.gb") if (m := re.fullmatch(r"batch_(\d+)\.gb", f.name))]
        self.next = max(nums, default=0) + 1

    def missing(self, accs: Iterable[str]) -> list[str]:
        return [a for a in accs if a not in self.where]

    def put_many(self, items: Iterable[tuple[str, str]]) -> int:
        recs = list(items)
        if not recs:
            return 0
        self.out.mkdir(parents=True, exist_ok=True)
        f = self.out / f"batch_{self.next:04d}.gb"
        self.next += 1
        tmp = f.with_suffix(".tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            for acc, text in recs:
                fh.write(text)
                self.where[acc] = f
        tmp.replace(f)
        return len(recs)

    def prune(self, keep: set[str]) -> int:
        """Remove records not in ``keep`` from the batch files; returns the number removed."""
        drop = {a for a in self.where if a not in keep}
        for f in sorted({self.where[a] for a in drop}):
            tmp = f.with_suffix(".tmp")
            n_left = 0
            with open(f, encoding="utf-8", errors="ignore") as src, open(tmp, "w", encoding="utf-8") as dst:
                for acc, text in iter_records(src):
                    if acc not in drop:
                        dst.write(text)
                        n_left += 1
            if n_left:
                tmp.replace(f)
            else:
                tmp.unlink()
                f.unlink()
        for a in drop:
            del self.where[a]
        return len(drop)

    def close(self) -> None:
        pass
