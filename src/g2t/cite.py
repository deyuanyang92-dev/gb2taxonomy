"""g2t-cite: find the published paper of GenBank records whose reference is still "Unpublished".

Many submitters never update GenBank after the paper appears, so the journal of a record stays
"Unpublished". For each distinct (title, authors) of such records, Crossref is searched by title and
first author; candidates are scored by title similarity (difflib ratio or word Jaccard, the larger)
and author-surname overlap:

    high    title_sim >= 0.90 and author_overlap >= 0.5
    medium  title_sim >= 0.75 and author_overlap >= 0.5 and >= 2 shared surnames (or all, if fewer)
    check   author_overlap >= 0.75, title_sim >= 0.5 and the taxon agrees (or title_sim >= 0.65):
            title changed before publication; or title_sim >= 0.90 with few shared authors
    low     anything else (not reported as a match)

A paper whose title does not name the record's taxon, when the GenBank title does, gets title_sim 0
(mitogenome notes share template titles: "The complete mitochondrial genome of ...").

Output (.xlsx): sheet "Curate" (high + medium rows, columns accession / Ref1Title / Ref1Authors /
Ref1Journal / curation_source) is table B for ``g2t-curate``; sheet "Candidates" lists every group
with scores and the Crossref record, "check" rows in yellow (copy them to "Curate" after checking).
A .csv output writes <stem>_curate.csv and <stem>_candidates.csv instead.

Crossref answers are cached in <output dir>/crossref_cache.json, so reruns query only new groups.
No e-mail is sent unless --mailto (or $CROSSREF_MAILTO) is given (Crossref "polite pool").

CLI:
    g2t-cite -i gb/mitogenome-header -o cite.xlsx            # GenBank files / folders
    g2t-cite -i final.csv -o cite.xlsx                       # a g2t table with Ref1..RefN columns
    g2t-curate -m matrix.xlsx -u cite.xlsx -o curated.xlsx   # apply the accepted matches
"""

from __future__ import annotations

import argparse
import html
import json
import os
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable, Iterator
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Optional

import pandas as pd

UA = "g2t-cite (gb2taxonomy; https://github.com/deyuanyang92-dev/gb2taxonomy)"
CROSSREF = "https://api.crossref.org/works"
# Crossref types that are papers (not "component" = figure/table of a paper, "peer-review", "dataset" ...)
PAPER_TYPES = {"journal-article", "proceedings-article", "posted-content", "book-chapter", "dissertation",
               "report", "book"}
TIERS = ["high", "medium", "check", "low", "none"]
CURATE_COLUMNS = ["accession", "Ref1Title", "Ref1Authors", "Ref1Journal", "curation_source"]
_NOT_PUBLISHED = ("Unpublished", "Submitted", "Direct Submission")

SearchFn = Callable[[str, str, Optional[int]], list[dict[str, Any]]]  # noqa: UP007 (runtime alias, Python 3.9)


# ------------------------------------------------------------------ text


def norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode()
    s = re.sub(r"<[^>]+>", " ", s)  # Crossref <i>…</i>
    s = re.sub(r"[^a-z0-9 ]", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def title_sim(a: str, b: str) -> float:
    a, b = norm(a), norm(b)
    if not a or not b:
        return 0.0
    ta, tb = set(a.split()), set(b.split())
    return round(max(SequenceMatcher(None, a, b).ratio(), len(ta & tb) / len(ta | tb)), 3)


def clean_title(s: str) -> str:
    """Crossref title -> plain text: tags removed, entities decoded (also double-escaped), spaces collapsed."""
    s = html.unescape(html.unescape(re.sub(r"<[^>]+>", " ", str(s or ""))))
    return re.sub(r"\s+", " ", s).strip()


_NOT_EPITHET = {"sp", "spp", "cf", "aff", "nr", "n", "gen", "nov", "var", "subsp", "x"}


def taxon_names(organisms: Iterable[str]) -> list[tuple[str, str]]:
    """(genus, epithet) of each organism ('Amynthas righii' -> ('amynthas', 'righii'); 'Megascolecidae sp.
    XX' -> ('megascolecidae', ''))."""
    out = []
    for o in organisms:
        w = norm(o).split()
        if not w:
            continue
        ep = w[1] if len(w) > 1 and w[1] not in _NOT_EPITHET else ""
        if (w[0], ep) not in out:
            out.append((w[0], ep))
    return out


def taxon_agrees(gb_title: str, cr_title: str, names: list[tuple[str, str]]) -> bool | None:
    """Does the paper's title name the taxon named in the GenBank title?

    True / False, or None when the GenBank title names none of the record's taxa. Mitogenome papers
    have template titles ("The complete mitochondrial genome of ..."), so title similarity alone pairs
    papers on different species; the taxon has to agree.
    """
    g, c = f" {norm(gb_title)} ", f" {norm(cr_title)} "
    named = [(ge, ep) for ge, ep in names if f" {ge} " in g]
    if not named:
        return None
    for ge, ep in named:
        if ep and f" {ep} " in g:
            if f" {ep} " in c and (f" {ge} " in c or f" {ge[0]} {ep} " in c):
                return True
        elif f" {ge} " in c:
            return True
    return False


def surnames(authors: str) -> list[str]:
    """GenBank AUTHORS ('Radashevsky,V.I., Pankova,V.V. and Lee,B.') -> normalised surnames."""
    out = []
    for p in re.split(r"\s+and\s+|\.,\s*", str(authors or "")):
        sn = norm(p.split(",")[0])
        if sn and sn not in out:
            out.append(sn)
    return out


def is_unpublished(journals: Iterable[str]) -> bool:
    js = [j for j in (str(x or "").strip() for x in journals) if j and j != "nan"]
    return bool(js) and not [j for j in js if not j.startswith(_NOT_PUBLISHED)]


def grade(title_s: float, overlap: float, shared: int | None = None, taxon_named: bool = False) -> str:
    """Tier of a candidate paper.

    shared = number of shared author surnames (None: not counted); one shared common surname
    (Wang, Li, Kim ...) is weak evidence, so medium needs two unless the record has one author.
    taxon_named = the GenBank title names the record's taxon and the paper's title names it too.
    """
    two = shared is None or shared >= 2 or overlap >= 0.99
    if title_s >= 0.9 and overlap >= 0.5:
        return "high"
    if title_s >= 0.75 and overlap >= 0.5 and two:
        return "medium"
    if overlap >= 0.75 and title_s >= 0.5 and (taxon_named or title_s >= 0.65):
        return "check"
    if title_s >= 0.9:
        return "check"
    return "low"


def genbank_journal(item: dict[str, Any]) -> str:
    """Crossref record -> GenBank JOURNAL style 'Zool Stud 65, e37 (2026)'."""
    short = clean_title((item.get("short-container-title") or [""])[0])
    full = clean_title((item.get("container-title") or [""])[0])
    j = full or short  # full title: unambiguous, and g2t-mitostats maps it to GenBank's abbreviation
    vol, page = item.get("volume", ""), item.get("page") or item.get("article-number", "")
    issue = item.get("issue", "")
    s = j
    if vol:
        s += f" {vol}" + (f" ({issue})" if issue else "")
    if page:
        s += f", {page}"
    y = cr_year(item)
    return f"{s} ({y})" if y else s


def cr_year(item: dict[str, Any]) -> int | None:
    for k in ("published-print", "published-online", "issued"):
        dp = (item.get(k) or {}).get("date-parts") or [[None]]
        if dp[0] and dp[0][0]:
            return int(dp[0][0])
    return None


def cr_authors(item: dict[str, Any]) -> str:
    """Crossref authors -> GenBank style 'Smith,A.B., Lee,C. and Wu,D.'."""
    names = []
    for a in item.get("author", []):
        fam = a.get("family") or a.get("name") or ""
        ini = "".join(f"{p[0]}." for p in re.split(r"[\s\-.]+", a.get("given", "")) if p)
        names.append(f"{fam},{ini}" if ini else fam)
    if len(names) > 1:
        return ", ".join(names[:-1]) + " and " + names[-1]
    return names[0] if names else ""


# ------------------------------------------------------------------ Crossref


def crossref_search(mailto: str = "") -> SearchFn:
    ua = UA + (f"; mailto:{mailto}" if mailto else "")

    def search(title: str, first_author: str, year: int | None) -> list[dict[str, Any]]:
        q = {"query.bibliographic": title, "rows": "5",
             "select": "DOI,title,author,container-title,short-container-title,issued,published-print,"
                       "published-online,type,volume,issue,page,article-number"}
        if first_author:
            q["query.author"] = first_author
        if year:
            q["filter"] = f"from-pub-date:{year - 1}"
        req = urllib.request.Request(CROSSREF + "?" + urllib.parse.urlencode(q), headers={"User-Agent": ua})
        err: Exception | None = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=30) as r:
                    return list(json.load(r)["message"]["items"])
            except Exception as e:  # noqa: BLE001
                err = e
                time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"Crossref: {err}")

    return search


def score_items(title: str, authors: str, items: list[dict[str, Any]],
                organisms: Iterable[str] = ()) -> dict[str, Any]:
    """Best Crossref record for a GenBank (title, authors): title similarity + 0.3 x author overlap.

    A paper whose title does not name the record's taxon (when the GenBank title does) gets
    title_sim 0. Figures/tables of a paper ("component") and other non-paper types are skipped.
    """
    sn = surnames(authors)
    names = taxon_names(organisms)
    best: dict[str, Any] = {}
    for it in items:
        if it.get("type") and it["type"] not in PAPER_TYPES:
            continue
        ct = clean_title((it.get("title") or [""])[0])
        agree = taxon_agrees(title, ct, names)
        ts = 0.0 if agree is False else title_sim(title, ct)
        cr_sn = {norm(a.get("family") or a.get("name") or "") for a in it.get("author", [])}
        shared = len(set(sn) & cr_sn)
        ov = round(shared / len(sn), 2) if sn else 0.0
        score = ts + 0.3 * ov
        if not best or score > best["_score"]:
            best = {"_score": score, "_shared": shared, "_agree": bool(agree),
                    "title_sim": ts, "author_overlap": ov, "doi": it.get("DOI", ""),
                    "cr_title": ct, "cr_authors": cr_authors(it), "cr_journal": genbank_journal(it),
                    "cr_year": cr_year(it), "cr_type": it.get("type", "")}
    if not best:
        return {"tier": "none"}
    best.pop("_score")
    best["tier"] = grade(best["title_sim"], best["author_overlap"], best.pop("_shared"), best.pop("_agree"))
    return best


def best_match(title: str, authors: str, year: int | None, search: SearchFn,
               organisms: Iterable[str] = ()) -> dict[str, Any]:
    sn = surnames(authors)
    return score_items(title, authors, search(title, sn[0] if sn else "", year), organisms)


# ------------------------------------------------------------------ input


def _from_genbank(paths: list[str]) -> Iterator[dict[str, Any]]:
    from g2t.mitostats import _gb_files, _records

    for rec in _records(_gb_files(paths)):
        refs = rec.annotations.get("references", [])
        yield {"accession": rec.id, "organism": rec.annotations.get("organism", ""),
               "refs": [(r.journal or "", r.title or "", r.authors or "") for r in refs]}


def _from_table(path: str) -> Iterator[dict[str, Any]]:
    df = pd.read_excel(path) if path.lower().endswith((".xlsx", ".xls")) else pd.read_csv(
        path, sep="\t" if path.lower().endswith(".tsv") else ",", low_memory=False, dtype=str)
    df = df.fillna("")
    acc = next((c for c in ("accession", "VERSION", "ACCESSION", "Accession") if c in df.columns), df.columns[0])
    org = next((c for c in ("organism", "ORGANISM", "Organism") if c in df.columns), None)
    n = max([int(m.group(1)) for c in df.columns if (m := re.match(r"Ref(\d+)Journal$", c))], default=0)
    if not n:
        raise SystemExit(f"{path}: no Ref1Journal .. RefNJournal columns")
    for _, r in df.iterrows():
        refs = [(str(r.get(f"Ref{i}Journal", "")), str(r.get(f"Ref{i}Title", "")), str(r.get(f"Ref{i}Authors", "")))
                for i in range(1, n + 1)]
        yield {"accession": str(r[acc]), "organism": str(r[org]) if org else "", "refs": refs}


def read_records(inputs: list[str]) -> Iterator[dict[str, Any]]:
    gb = [p for p in inputs if not p.lower().endswith((".csv", ".tsv", ".xlsx", ".xls"))]
    if gb:
        yield from _from_genbank(gb)
    for p in inputs:
        if p not in gb:
            yield from _from_table(p)


def unpublished_groups(records: Iterable[dict[str, Any]]) -> pd.DataFrame:
    """Records whose references are all Unpublished/Submitted, grouped by (title, authors)."""
    rows = []
    for r in records:
        refs = r["refs"]
        if not is_unpublished(j for j, _, _ in refs):
            continue
        title, authors = next(((t, a) for j, t, a in refs
                               if j.startswith("Unpublished") and t and t != "Direct Submission"), ("", ""))
        if not title:
            continue
        years = [int(m.group(1)) for j, _, _ in refs if (m := re.search(r"Submitted \(\d{2}-\w{3}-(\d{4})\)", j))]
        rows.append({"accession": r["accession"], "organism": r["organism"], "title": title, "authors": authors,
                     "submit_year": max(years) if years else None})
    cols = ["title", "authors", "submit_year", "n_records", "accessions", "organisms"]
    if not rows:
        return pd.DataFrame(columns=cols)
    d = pd.DataFrame(rows)
    g = d.groupby(["title", "authors"], sort=False).agg(
        submit_year=("submit_year", "max"), n_records=("accession", "size"),
        accessions=("accession", "; ".join),
        organisms=("organism", lambda s: "; ".join(dict.fromkeys(s))),
    ).reset_index()
    return g[cols]


# ------------------------------------------------------------------ run


def _key(title: str, authors: str) -> str:
    return norm(title) + "|" + "/".join(surnames(authors)[:3])


CACHE_VERSION = 2


def match_groups(groups: pd.DataFrame, search: SearchFn, cache_path: Path | None = None,
                 delay: float = 0.2, progress: bool = True) -> pd.DataFrame:
    """Score each (title, authors) group. The raw Crossref answers are cached, so a rerun (or a change
    of the scoring rules) queries only groups not seen before."""
    cache: dict[str, list[dict[str, Any]]] = {}
    if cache_path and cache_path.exists():
        saved = json.loads(cache_path.read_text(encoding="utf-8"))
        if isinstance(saved, dict) and saved.get("version") == CACHE_VERSION:
            cache = saved["items"]

    def save() -> None:
        if cache_path:
            cache_path.write_text(json.dumps({"version": CACHE_VERSION, "items": cache}, ensure_ascii=False),
                                  encoding="utf-8")

    out = []
    todo = sum(_key(t, a) not in cache for t, a in zip(groups["title"], groups["authors"]))
    done = 0
    for _, g in groups.iterrows():
        k = _key(g["title"], g["authors"])
        if k not in cache:
            y = int(g["submit_year"]) if pd.notna(g["submit_year"]) else None
            sn = surnames(g["authors"])
            try:
                cache[k] = search(g["title"], sn[0] if sn else "", y)
            except Exception as e:  # noqa: BLE001  (network: report, keep going, retry next run)
                print(f"   {g['title'][:60]}: {e}", file=sys.stderr)
                out.append({**g.to_dict(), "tier": "error"})
                continue
            done += 1
            if done % 20 == 0:
                save()
                if progress:
                    print(f"   Crossref {done}/{todo}", file=sys.stderr)
            time.sleep(delay)
        out.append({**g.to_dict(), **score_items(g["title"], g["authors"], cache[k],
                                                  str(g.get("organisms", "")).split("; "))})
    save()
    t = pd.DataFrame(out)
    t["tier"] = pd.Categorical(t["tier"], TIERS + ["error"], ordered=True) if len(t) else []
    return t.sort_values(["tier", "n_records"], ascending=[True, False]).reset_index(drop=True) if len(t) else t


def curate_table(cand: pd.DataFrame, tiers: tuple[str, ...] = ("high", "medium")) -> pd.DataFrame:
    """One row per accession for g2t-curate (table B)."""
    rows = []
    for _, c in cand[cand["tier"].isin(tiers)].iterrows():
        src = f"Crossref {c['tier']} match doi:{c['doi']} (title_sim {c['title_sim']}, authors {c['author_overlap']})"
        for a in str(c["accessions"]).split("; "):
            rows.append([a, c["cr_title"], c["cr_authors"], c["cr_journal"], src])
    return pd.DataFrame(rows, columns=CURATE_COLUMNS)


def write_output(cand: pd.DataFrame, path: Path) -> None:
    cur = curate_table(cand)
    cols = ["tier", "title_sim", "author_overlap", "n_records", "title", "cr_title", "authors", "cr_authors",
            "cr_journal", "cr_year", "doi", "cr_type", "submit_year", "organisms", "accessions"]
    cand = cand[[c for c in cols if c in cand.columns]].copy()
    if "doi" in cand:
        cand["doi"] = [f"https://doi.org/{d}" if isinstance(d, str) and d else "" for d in cand["doi"]]
    if path.suffix.lower() != ".xlsx":
        cur.to_csv(path.with_name(path.stem + "_curate.csv"), index=False)
        cand.to_csv(path.with_name(path.stem + "_candidates.csv"), index=False)
        return
    from openpyxl.styles import Font, PatternFill

    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        cur.to_excel(xw, sheet_name="Curate", index=False)
        cand.to_excel(xw, sheet_name="Candidates", index=False)
        yellow = PatternFill("solid", fgColor="FFF2CC")
        ws = xw.sheets["Candidates"]
        for i, tier in enumerate(cand["tier"].astype(str), start=2):
            if tier == "check":
                for cell in ws[i]:
                    cell.fill = yellow
        for w in xw.sheets.values():
            w.freeze_panes = "A2"
            for cell in w[1]:
                cell.font = Font(bold=True)
            for col in w.columns:
                w.column_dimensions[col[0].column_letter].width = min(
                    60, max(10, *(len(str(c.value or "")) for c in col[:50])))


def run(inputs: list[str], output: str, search: SearchFn | None = None, delay: float = 0.2,
        progress: bool = True) -> pd.DataFrame:
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    groups = unpublished_groups(read_records(inputs))
    if progress:
        print(f"   {int(groups['n_records'].sum()) if len(groups) else 0} unpublished records, "
              f"{len(groups)} distinct (title, authors)", file=sys.stderr)
    if not len(groups):
        cand = pd.DataFrame(columns=["tier", "accessions"])
    else:
        cand = match_groups(groups, search or crossref_search(os.environ.get("CROSSREF_MAILTO", "")),
                            cache_path=out.parent / "crossref_cache.json", delay=delay, progress=progress)
    write_output(cand, out)
    return cand


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="g2t-cite",
        description="Find the published paper (Crossref) of GenBank records still marked 'Unpublished'.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
        "  g2t-cite -i gb/mitogenome-header -o cite.xlsx\n"
        "  g2t-cite -i final.csv -o cite.xlsx\n"
        "  g2t-curate -m matrix.xlsx -u cite.xlsx -o curated.xlsx   # apply sheet 'Curate'\n"
        "  g2t-mitostats -i gb/mitogenome-header -o mito --citations cite.xlsx",
    )
    p.add_argument("-i", "--input", nargs="+", required=True,
                   help="GenBank files/folders, or g2t tables (.csv/.tsv/.xlsx) with Ref1Journal.. columns")
    p.add_argument("-o", "--output", required=True, help="cite.xlsx (or .csv -> _curate.csv + _candidates.csv)")
    p.add_argument("--mailto", default=os.environ.get("CROSSREF_MAILTO", ""),
                   help="your e-mail for the Crossref polite pool (default: $CROSSREF_MAILTO; none sent otherwise)")
    p.add_argument("--delay", type=float, default=0.2, help="seconds between Crossref queries (default 0.2)")
    p.add_argument("-q", "--quiet", action="store_true")
    a = p.parse_args(argv)
    cand = run(a.input, a.output, search=crossref_search(a.mailto), delay=a.delay, progress=not a.quiet)
    counts = cand["tier"].astype(str).value_counts().to_dict() if len(cand) else {}
    print(f"g2t-cite: {a.output}  " + ", ".join(f"{t} {counts[t]}" for t in TIERS + ["error"] if t in counts))
    return 0


if __name__ == "__main__":
    sys.exit(main())
