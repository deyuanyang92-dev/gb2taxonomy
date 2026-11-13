"""g2t-mitostats: sequencing technology, assembly programs and publication of mitochondrial genomes.

Reads GenBank files (full records or ``g2t-download --header-only`` headers), takes from each record
the "##Assembly-Data##" comment (Sequencing Technology, Assembly Method), the references and the
taxon, and writes

* ``mitostats.csv``            one row per record, raw GenBank text next to every derived value;
* ``counts_by_year.csv``       records per year and sequencing category (and cumulative);
* ``counts_by_assembler.csv``  records per assembler group and sequencing category;
* ``counts_by_journal.csv``    published records per journal;
* ``fig_seqtech_by_year.{png,pdf}`` cumulative area chart + pie (cf. Yang, PhD thesis, Fig. 2-12);
* ``fig_alluvial.{png,pdf}``   year -> sequencing category -> assembler (cf. Fig. 2-13);
* ``mitostats_report.html``    one offline, interactive page (filters by year, category, assembler, taxon).

Classification is conservative: a record without an Assembly-Data comment is "Unknown", never assumed
to be Sanger. Rules and spellings are in ``data/mito_methods.yaml``.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
import time
from collections.abc import Iterable, Iterator
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

CATEGORIES = ["Sanger", "NGS", "Long-read", "Unknown"]
RANKS = ["phylum", "class", "order", "family", "genus"]
NOT_STATED = "Not stated"
MULTI = "Multi-assembly"
OTHER = "Other"


# ------------------------------------------------------------------ rules


@lru_cache(maxsize=1)
def _rules() -> dict[str, Any]:
    import yaml

    text = resources.files("g2t").joinpath("data/mito_methods.yaml").read_text(encoding="utf-8")
    raw = yaml.safe_load(text)

    def comp(pats: list[str]) -> re.Pattern[str]:
        return re.compile("|".join(f"(?:{p})" for p in pats), re.IGNORECASE)

    return {
        "platforms": [(p["name"], p["category"], comp(p["patterns"])) for p in raw["platforms"]],
        "assemblers": [(a["name"], comp(a["patterns"])) for a in raw["assemblers"]],
        "not_assembly": comp(raw["not_assembly"]),
    }


def sequencing_platforms(raw: str) -> list[str]:
    """Platforms named in a Sequencing Technology value, in rule order."""
    raw = str(raw or "")
    return [name for name, _, rx in _rules()["platforms"] if rx.search(raw)]


def sequencing_category(raw: str) -> str:
    """Sanger / NGS / Long-read / Unknown. Long-read wins (hybrid assemblies), then NGS, then Sanger."""
    raw = str(raw or "")
    cats = {cat for _, cat, rx in _rules()["platforms"] if rx.search(raw)}
    for c in ("Long-read", "NGS", "Sanger"):
        if c in cats:
            return c
    return "Unknown"


_SPLIT = re.compile(r"\s*;\s*|\s*,\s*(?=[A-Za-z])|\s+and\s+|\s*\+\s*|\s*&\s*", re.IGNORECASE)


def assemblers(raw: str) -> tuple[list[str], list[str]]:
    """(recognised assembly programs, unrecognised entries) of an Assembly Method value.

    Quality-control, annotation and mapping tools are dropped (``not_assembly`` in the rules).
    """
    rules = _rules()
    names: list[str] = []
    other: list[str] = []
    for item in _SPLIT.split(str(raw or "")):
        item = item.strip()
        if not item:
            continue
        hit = [n for n, rx in rules["assemblers"] if rx.search(item)]
        if hit:
            names += [n for n in hit if n not in names]
        elif (
            not rules["not_assembly"].search(item)
            and re.search(r"[A-Za-z0-9]", item)
            and not re.fullmatch(r"v\.?\s*[\w.]+", item, re.IGNORECASE)  # a stray version ("v. 4.3")
        ):
            other.append(item)
    return names, other


def assembler_group(names: list[str], other: list[str]) -> str:
    """One label per record: the program, "Multi-assembly", "Other" or "Not stated"."""
    if len(names) > 1:
        return MULTI
    if names:
        return names[0]
    return OTHER if other else NOT_STATED


# ------------------------------------------------------------------ references

_NOT_PUBLISHED = re.compile(r"^\s*(unpublished|submitted|direct submission)", re.IGNORECASE)
_SUBMITTED = re.compile(r"Submitted \(\d{1,2}-[A-Z]{3}-(\d{4})\)", re.IGNORECASE)
_YEAR = re.compile(r"\((\d{4})\)")


def journal_name_year(journal: str) -> tuple[str, int | None]:
    """'Zool Stud 65, e37 (2026)' -> ('Zool Stud', 2026); unpublished/submitted -> ('', None)."""
    j = str(journal or "").strip()
    if not j or _NOT_PUBLISHED.match(j):
        return "", None
    years = _YEAR.findall(j)
    year = int(years[-1]) if years else None
    name = re.split(r"\s+\d|\s*\(|,", j, maxsplit=1)[0]
    # GenBank writes one journal several ways ("Mol. Phylogenet. Evol." / "Mol Phylogenet Evol"): no dots
    name = re.sub(r"\s+", " ", name.replace(".", " ")).strip()
    return name, year


def publication(refs: list[tuple[str, str]], create_date: str = "") -> dict[str, Any]:
    """Journal, status and year from a record's references [(JOURNAL, TITLE), ...].

    year = publication year, else the year of the direct submission, else of the record's creation.
    """
    journal, pub_year = "", None
    submit_year = None
    for j, _title in refs:
        m = _SUBMITTED.search(str(j or ""))
        if m and submit_year is None:
            submit_year = int(m.group(1))
        name, y = journal_name_year(j)
        if name and not journal:
            journal, pub_year = name, y
    status = "published" if journal else "unpublished"
    if journal.lower().startswith("thesis"):
        status = "thesis"
    created = int(create_date[:4]) if re.match(r"\d{4}", str(create_date or "")) else None
    if pub_year:
        year, source = pub_year, "publication"
    elif submit_year:
        year, source = submit_year, "submission"
    elif created:
        year, source = created, "record created"
    else:
        year, source = None, ""
    return {"journal": journal, "pub_status": status, "pub_year": pub_year, "submit_year": submit_year,
            "year": year, "year_source": source}


# ------------------------------------------------------------------ taxonomy


def fetch_lineage(taxids: Iterable[str], cache_path: Path | None = None) -> dict[str, dict[str, str]]:
    """Phylum..genus for NCBI TaxIDs (Entrez taxonomy, 200 per request), cached as JSON."""
    cache: dict[str, dict[str, str]] = {}
    if cache_path and cache_path.exists():
        cache = json.loads(cache_path.read_text(encoding="utf-8"))
    todo = sorted({t for t in taxids if t and t not in cache})
    if todo:
        from Bio import Entrez

        from g2t.download import _call, setup_entrez

        delay = setup_entrez()
        for i in range(0, len(todo), 200):
            part = todo[i : i + 200]
            try:
                recs = _call(Entrez.efetch, True, db="taxonomy", id=",".join(part), retmode="xml")
            except Exception as e:  # noqa: BLE001 - taxonomy is optional; columns stay empty
                logger.warning(f"NCBI taxonomy lookup failed for {len(part)} TaxIDs: {e}")
                continue
            for rec in recs:
                nodes = list(rec.get("LineageEx", [])) + [rec]
                lin = {str(n.get("Rank")): str(n.get("ScientificName")) for n in nodes if n.get("Rank") in RANKS}
                for t in [str(rec.get("TaxId", ""))] + [str(x) for x in rec.get("AkaTaxIds", [])]:
                    cache[t] = lin
            time.sleep(delay)
        if cache_path:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
    return cache


# ------------------------------------------------------------------ records -> table


def _gb_files(paths: Iterable[str | Path]) -> list[Path]:
    out: list[Path] = []
    for p in map(Path, paths):
        if p.is_dir():
            out += sorted(f for f in p.rglob("*") if re.search(r"\.(gb|gbk|gbff|genbank)(\.gz)?$", f.name))
        elif p.exists():
            out.append(p)
    return out


def _summaries(files: list[Path]) -> dict[str, dict[str, str]]:
    """summary.tsv files written by g2t-download --header-only (true length and dates)."""
    out: dict[str, dict[str, str]] = {}
    for d in {f.parent for f in files}:
        s = d / "summary.tsv"
        if s.exists():
            with open(s, encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    out[row["accession"]] = row
    return out


def _records(files: list[Path]) -> Iterator[Any]:
    import warnings

    from Bio import BiopythonParserWarning, SeqIO

    from g2t._extract.files import open_genbank

    for f in files:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open_genbank(str(f)) as fh:
                yield from SeqIO.parse(fh, "genbank")


def _record_row(rec: Any, summary: dict[str, str]) -> dict[str, Any]:
    from g2t._extract.parse import extract_assembly_from_record

    asm, _src = extract_assembly_from_record(rec)
    asm = asm or {}
    taxid = ""
    for f in rec.features:
        if f.type == "source":
            for x in f.qualifiers.get("db_xref", []):
                if x.startswith("taxon:"):
                    taxid = x.split(":", 1)[1]
            break
    taxid = taxid or summary.get("taxid", "")
    refs = [(r.journal or "", r.title or "") for r in rec.annotations.get("references", [])]
    pub = publication(refs, summary.get("create_date", ""))
    first = next(((j, t, a) for j, t, a in ((r.journal, r.title, r.authors) for r in
                  rec.annotations.get("references", [])) if t and t != "Direct Submission"), ("", "", ""))
    st_raw = asm.get("Sequencing Technology", "")
    am_raw = asm.get("Assembly Method", "")
    names, other = assemblers(am_raw)
    length = summary.get("length") or (len(rec.seq) if len(rec.seq) > 1 else None)
    return {
        "accession": rec.id,
        "organism": rec.annotations.get("organism", ""),
        "taxid": taxid,
        "length": int(length) if length else None,
        "create_date": summary.get("create_date", ""),
        "lineage": "; ".join(rec.annotations.get("taxonomy", [])),
        "seq_technology_raw": st_raw,
        "seq_platforms": "; ".join(sequencing_platforms(st_raw)),
        "seq_category": sequencing_category(st_raw),
        "assembly_method_raw": am_raw,
        "assemblers": "; ".join(names),
        "assembler_group": assembler_group(names, other),
        "assembly_other": "; ".join(other),
        "genome_coverage": asm.get("Genome Coverage", ""),
        "has_assembly_data": bool(asm),
        **pub,
        "journal_raw": first[0] or "",
        "title": first[1] or "",
        "authors": first[2] or "",
    }


def build_table(paths: Iterable[str | Path], lineage: dict[str, dict[str, str]] | None = None,
                cache_path: Path | None = None) -> pd.DataFrame:
    files = _gb_files(paths)
    if not files:
        raise FileNotFoundError(f"no GenBank files in {list(paths)}")
    summ = _summaries(files)
    seen: set[str] = set()
    rows = []
    for rec in _records(files):
        if rec.id in seen:
            continue
        seen.add(rec.id)
        rows.append(_record_row(rec, summ.get(rec.id, {})))
    t = pd.DataFrame(rows)
    lin = lineage if lineage is not None else fetch_lineage(t["taxid"].astype(str), cache_path)
    for r in RANKS:
        t.insert(t.columns.get_loc("lineage"), r, [lin.get(str(x), {}).get(r, "") for x in t["taxid"]])
    t["year"] = t["year"].astype("Int64")
    for c in ("pub_year", "submit_year", "length"):
        t[c] = t[c].astype("Int64")
    return t


# ------------------------------------------------------------------ summaries


def counts_by_year(t: pd.DataFrame) -> pd.DataFrame:
    c = t.dropna(subset=["year"]).pivot_table(index="year", columns="seq_category", values="accession",
                                              aggfunc="count", fill_value=0)
    c = c.reindex(columns=CATEGORIES, fill_value=0)
    if len(c):
        c = c.reindex(range(int(c.index.min()), int(c.index.max()) + 1), fill_value=0)
    c.index.name = "year"
    for k in CATEGORIES:
        c[f"{k}_cumulative"] = c[k].cumsum()
    return c.reset_index()


def counts_by_assembler(t: pd.DataFrame) -> pd.DataFrame:
    c = t.pivot_table(index="assembler_group", columns="seq_category", values="accession", aggfunc="count",
                      fill_value=0).reindex(columns=CATEGORIES, fill_value=0)
    c["total"] = c.sum(axis=1)
    return c.sort_values("total", ascending=False).reset_index()


def counts_by_journal(t: pd.DataFrame) -> pd.DataFrame:
    p = t[t["pub_status"] == "published"]
    return (p.groupby("journal").size().rename("records").sort_values(ascending=False).reset_index())


# ------------------------------------------------------------------ figures

COLORS = {"Sanger": "#2a78d6", "NGS": "#eb6834", "Long-read": "#1baf7a", "Unknown": "#a3a29b"}


def fig_seqtech_by_year(cy: pd.DataFrame, out: Path, include_unknown: bool = True) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cats = [c for c in CATEGORIES if include_unknown or c != "Unknown"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ys = [cy[f"{c}_cumulative"] for c in cats]
    ax.stackplot(cy["year"], *ys, labels=cats, colors=[COLORS[c] for c in cats], alpha=0.9,
                 edgecolor="white", linewidth=0.8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Mitochondrial genomes (cumulative)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e5e4df", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, bbox_to_anchor=(0.0, 1.0))
    tot = [int(cy[c].sum()) for c in cats]
    ins = ax.inset_axes((0.30, 0.42, 0.32, 0.55))
    ins.pie(tot, colors=[COLORS[c] for c in cats], startangle=90, counterclock=False,
            wedgeprops={"edgecolor": "white", "linewidth": 1.5})
    ins.set_title(f"n = {sum(tot)}", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_seqtech_by_year.{ext}", dpi=300)
    plt.close(fig)


def _top(series: pd.Series, n: int, rest: str = OTHER) -> pd.Series:
    keep = series.value_counts().index[:n]
    return series.where(series.isin(keep), rest)


def fig_alluvial(t: pd.DataFrame, out: Path, top_assemblers: int = 15) -> None:
    """Year -> sequencing category -> assembler group, bands coloured by sequencing category."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MPath

    d = t.dropna(subset=["year"]).copy()
    d["year"] = d["year"].astype(int).astype(str)
    d["asm"] = _top(d["assembler_group"], top_assemblers)
    cols = ["year", "seq_category", "asm"]
    orders = [
        sorted(d["year"].unique()),
        [c for c in CATEGORIES if c in set(d["seq_category"])],
        [a for a in d["asm"].value_counts().index if a not in (OTHER, NOT_STATED)]
        + [a for a in (OTHER, NOT_STATED) if a in set(d["asm"])],
    ]
    n = len(d)
    gap = n * 0.004
    x = [0.0, 1.0, 2.0]
    w = 0.12
    fig, ax = plt.subplots(figsize=(9, 9))
    # node positions (top -> bottom)
    pos: list[dict[str, list[float]]] = []
    for ci, col in enumerate(cols):
        sizes = d[col].value_counts()
        y = 0.0
        p: dict[str, list[float]] = {}
        for k in orders[ci]:
            p[k] = [y, y + sizes[k]]
            y += sizes[k] + gap
        pos.append(p)
        for k, (y0, y1) in p.items():
            ax.add_patch(plt.Rectangle((x[ci] - w / 2, -y1), w, y1 - y0, facecolor="#3d3d3a", edgecolor="white"))
            if y1 - y0 > n * 0.008:
                ha, dx = ("right", -w / 2 - 0.02) if ci == 0 else ("left", w / 2 + 0.02)
                ax.text(x[ci] + dx, -(y0 + y1) / 2, k, ha=ha, va="center", fontsize=7)
    # bands
    for ci in range(2):
        a, b = cols[ci], cols[ci + 1]
        cur_a = {k: v[0] for k, v in pos[ci].items()}
        cur_b = {k: v[0] for k, v in pos[ci + 1].items()}
        flows: dict[tuple[str, str, str], int] = {}
        for ka, kb, cat in zip(d[a], d[b], d["seq_category"]):
            flows[(ka, kb, cat)] = flows.get((ka, kb, cat), 0) + 1
        for ka in orders[ci]:
            for kb in orders[ci + 1]:
                for cat in CATEGORIES:
                    v = flows.get((ka, kb, cat), 0)
                    if not v:
                        continue
                    ya, yb = cur_a[ka], cur_b[kb]
                    x0, x1 = x[ci] + w / 2, x[ci + 1] - w / 2
                    xm = (x0 + x1) / 2
                    verts = [(x0, -ya), (xm, -ya), (xm, -yb), (x1, -yb), (x1, -(yb + v)), (xm, -(yb + v)),
                             (xm, -(ya + v)), (x0, -(ya + v)), (x0, -ya)]
                    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
                             MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
                    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=COLORS[cat], alpha=0.55, lw=0))
                    cur_a[ka] += v
                    cur_b[kb] += v
    total_h = max(max(v[1] for v in p.values()) for p in pos)
    ax.set_xlim(-0.6, 2.7)
    ax.set_ylim(-total_h * 1.01, total_h * 0.03)
    for xi, lab in zip(x, ["Year", "Sequencing", "Assembly"]):
        ax.text(xi, total_h * 0.015, lab, ha="center", va="bottom", fontsize=9, weight="bold")
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_alluvial.{ext}", dpi=300)
    plt.close(fig)


# ------------------------------------------------------------------ HTML report

REPORT_COLUMNS = ["accession", "organism", "phylum", "class", "order", "family", "year", "year_source",
                  "seq_category", "seq_platforms", "assembler_group", "assemblers", "journal", "pub_status",
                  "journal_source", "length", "seq_technology_raw", "assembly_method_raw", "title"]


def _plotly_js() -> str:
    """plotly.min.js of the installed plotly package (inlined: the report works offline)."""
    import plotly

    return (Path(plotly.__file__).parent / "package_data" / "plotly.min.js").read_text(encoding="utf-8")


def write_report(t: pd.DataFrame, path: Path, title: str = "") -> None:
    """One self-contained offline HTML page: data embedded as JSON, plotly.js inlined (~5 MB)."""
    tpl = resources.files("g2t").joinpath("data/mitostats_report.html").read_text(encoding="utf-8")
    recs = t[[c for c in REPORT_COLUMNS if c in t.columns]].astype(object).where(t.notna(), None)
    data = {"title": title or "Mitochondrial genomes", "generated": time.strftime("%Y-%m-%d"),
            "categories": CATEGORIES, "colors": COLORS, "records": recs.to_dict(orient="records")}
    plotly_tag = f"<script>{_plotly_js()}</script>"
    payload = json.dumps(data, ensure_ascii=False, default=str).replace("</", "<\\/")
    html = tpl.replace("__G2T_TITLE__", data["title"]).replace("__G2T_PLOTLY__", plotly_tag)
    html = html.replace("__G2T_DATA__", payload)
    path.write_text(html, encoding="utf-8")


# ------------------------------------------------------------------ run / CLI


_JOURNAL_STOP = {"part", "of", "the", "and", "&", "section", "series"}


def _jwords(name: str) -> list[str]:
    return [w for w in re.sub(r"[^a-z0-9 ]", " ", name.lower()).split() if w not in _JOURNAL_STOP]


def _prefix_subsequence(short: list[str], long: list[str]) -> bool:
    i = 0
    for w in short:
        while i < len(long) and not long[i].startswith(w):
            i += 1
        if i == len(long):
            return False
        i += 1
    return True


def same_journal(a: str, b: str) -> bool:
    """'Mitochondrial DNA B Resour' ~ 'Mitochondrial DNA Part B' ~ 'Mitochondrial DNA Part B: Resources';
    'Mol Phylogenet Evol' ~ 'Molecular Phylogenetics and Evolution'. Word-prefix match in order, either way."""
    wa, wb = _jwords(a), _jwords(b)
    if not wa or not wb or not (wa[0].startswith(wb[0]) or wb[0].startswith(wa[0])):
        return False
    if len(wa) == 1 or len(wb) == 1:  # "Gene" is not "Genes", "Diversity" not "Diversity and Distributions"
        return wa == wb
    return _prefix_subsequence(wa, wb) or _prefix_subsequence(wb, wa)


def apply_citations(t: pd.DataFrame, path: str | Path) -> pd.DataFrame:
    """Journal and publication year of unpublished records from g2t-cite (sheet 'Curate' / *_curate.csv).

    Only records still unpublished in GenBank change; journal_source says where the journal came from.
    """
    p = Path(path)
    if p.suffix.lower() == ".xlsx":
        c = pd.read_excel(p, sheet_name="Curate", dtype=str)
    else:
        c = pd.read_csv(p.with_name(p.stem + "_curate.csv") if not p.stem.endswith("_curate") else p, dtype=str)
    found = {str(r["accession"]): str(r["Ref1Journal"]) for _, r in c.fillna("").iterrows()}
    t = t.copy()
    t["journal_source"] = t["pub_status"].map(lambda s: "GenBank" if s != "unpublished" else "")
    # spell a Crossref journal the way GenBank records already spell it (commonest form first)
    known = list(t.loc[t["pub_status"] == "published", "journal"].value_counts().index)
    n = 0
    for i, r in t.iterrows():
        j = found.get(str(r["accession"])) or found.get(str(r["accession"]).split(".")[0])
        if r["pub_status"] != "unpublished" or not j:
            continue
        name, y = journal_name_year(j)
        if not name:
            continue
        match = [k for k in known if same_journal(k, name)]  # most specific GenBank spelling wins
        name = max(match, key=lambda k: len(_jwords(k))) if match else name
        t.at[i, "journal"], t.at[i, "pub_status"], t.at[i, "journal_source"] = name, "published", "Crossref"
        if y:
            t.at[i, "pub_year"], t.at[i, "year"], t.at[i, "year_source"] = y, y, "publication (Crossref)"
        n += 1
    logger.info("citations: %d unpublished records given a journal from %s", n, p.name)
    return t


def run(inputs: list[str], out_dir: str, lineage: dict[str, dict[str, str]] | None = None,
        figures: bool = True, title: str = "", top_assemblers: int = 15,
        citations: str | None = None) -> pd.DataFrame:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t = build_table(inputs, lineage=lineage, cache_path=out / "taxonomy_lineage.json")
    if citations:
        t = apply_citations(t, citations)
    t.to_csv(out / "mitostats.csv", index=False)
    cy = counts_by_year(t)
    cy.to_csv(out / "counts_by_year.csv", index=False)
    counts_by_assembler(t).to_csv(out / "counts_by_assembler.csv", index=False)
    counts_by_journal(t).to_csv(out / "counts_by_journal.csv", index=False)
    if figures:
        try:
            fig_seqtech_by_year(cy, out)
            fig_alluvial(t, out, top_assemblers)
        except ImportError:
            logger.warning("matplotlib not installed: figures skipped (pip install matplotlib)")
    try:
        write_report(t, out / "mitostats_report.html", title)
    except ImportError:
        logger.warning("plotly not installed: HTML report skipped (pip install plotly)")
    return t


def describe(t: pd.DataFrame) -> list[str]:
    n = len(t)
    lines = [f"{n} records; {int(t['has_assembly_data'].sum())} with an Assembly-Data comment"]
    lines.append("sequencing: " + ", ".join(f"{k} {int((t['seq_category'] == k).sum())}" for k in CATEGORIES))
    top = t["assembler_group"].value_counts().head(8)
    lines.append("assembly:   " + ", ".join(f"{k} {v}" for k, v in top.items()))
    pub = int((t["pub_status"] == "published").sum())
    nj = t.loc[t["pub_status"] == "published", "journal"].nunique()
    lines.append(f"published:  {pub} ({100 * pub / max(n, 1):.0f}%), {nj} journals")
    if t["year"].notna().any():
        lines.append(f"years:      {int(t['year'].min())}-{int(t['year'].max())}")
    return lines


EPILOG = """examples:
  g2t-download -t Annelida -o gb --mitogenome --header-only
  g2t-mitostats -i gb/mitogenome-header -o annelida_mito
  open annelida_mito/mitostats_report.html
  g2t-cite -i gb/mitogenome-header -o annelida_mito/cite.xlsx      # papers of "Unpublished" records
  g2t-mitostats -i gb/mitogenome-header -o annelida_mito --citations annelida_mito/cite.xlsx

Records without an Assembly-Data comment are counted as "Unknown" (not as Sanger).
Spellings and rules: g2t/data/mito_methods.yaml."""


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="g2t-mitostats",
        description="Sequencing technology, assembly programs and publications of mitochondrial genomes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EPILOG,
    )
    p.add_argument("-i", "--input", nargs="+", required=True, metavar="PATH",
                   help="GenBank files or directories (e.g. gb/mitogenome-header)")
    p.add_argument("-o", "--output", required=True, metavar="DIR", help="output directory")
    p.add_argument("--title", default="", help="title of the HTML report")
    p.add_argument("--top", type=int, default=15, metavar="N", help="assemblers shown in the alluvial figure")
    p.add_argument("--no-figures", action="store_true", help="skip the PNG/PDF figures")
    p.add_argument("--no-taxonomy", action="store_true", help="do not look up phylum..genus on NCBI")
    p.add_argument("--citations", metavar="FILE",
                   help="g2t-cite output (cite.xlsx): journals of records still 'Unpublished' in GenBank")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    t = run(a.input, a.output, lineage={} if a.no_taxonomy else None, figures=not a.no_figures,
            title=a.title, top_assemblers=a.top, citations=a.citations)
    for line in describe(t):
        print(line)
    print(f"output: {a.output}/mitostats.csv, mitostats_report.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())
