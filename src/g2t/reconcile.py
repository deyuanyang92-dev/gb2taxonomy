#!/usr/bin/env python3
"""
Voucher Reconciler (Step 3b)
============================
Merges specimen groups whose vouchers are written differently but refer to the
same specimen, e.g. ``COI_ZMMU_MSU_WS399`` / ``28S_ZMMU_WS399`` / ``WS399``.

Exact string matching (Step 3) splits such specimens into one row per gene.
This step decides merges from evidence in the GenBank records themselves,
never from the voucher string alone:

1. Candidates: groups of the *same organism* whose vouchers share a core
   identifier (gene names and collection acronyms stripped; >= 4 characters).
   The same core in different organisms is reported, never merged.
2. Evidence for each candidate pair:
     strong    - a shared publication title (not "Direct Submission");
                 identical full collection date; coordinates within 0.01 deg
     moderate  - same first author; same collector; same detailed locality;
                 coordinates within 0.5 deg
     conflict  - incompatible collection dates; different countries;
                 coordinates > 0.5 deg apart
3. Decision: any conflict -> not merged; >= 1 strong -> merge ("high");
   >= 2 moderate -> merge ("medium"); otherwise not merged. If both groups
   already contain the same gene, only strong evidence can merge them.
   Merges are applied greedily (strongest first) and a merge is skipped if it
   would put two conflicting groups in one specimen.

Outputs (in ``output_dir``):
  reconciled_species_voucher.csv  input rows + species_voucher_g2t, voucher_core,
                                  match_basis (exact/reconciled),
                                  match_confidence, match_evidence;
                                  species_voucher_new holds the reconciled key
  reconcile_report.csv            one row per candidate pair with the decision
"""

from __future__ import annotations

import argparse
import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from g2t.utils import StepResult, read_table

GENE_TOKENS = (r"COI|CO1|COX1|COII|COX2|COIII|COX3|CYTB|COB|ND\d|NAD\d|ATP\d|12S|16S|18S|28S|5\.8S|"
               r"H3|ITS1?|ITS2|EF1A?|EF1ALPHA|RRNL|RRNS|LSU|SSU|MTDNA|MITO")
GENE_PREFIX = re.compile(rf"^(?:{GENE_TOKENS})[\s_\-:]+", re.I)
GENE_SUFFIX = re.compile(rf"[\s_\-:]+(?:{GENE_TOKENS})$", re.I)
CORE = re.compile(r"(?:^|[^A-Za-z0-9])([A-Za-z]{0,6})[\s_\-]?((?:(?:1[89]|20)\d{2}[\s_.\-])?\d{2,8}[A-Za-z]?)$")
ID_TOKEN = re.compile(r"[A-Za-z]{1,6}\d{3,8}[A-Za-z]?")
GENE_FULL = re.compile(rf"(?:{GENE_TOKENS})", re.I)
TRAILING = re.compile(r"(\s*\([^()]*\)|[\s.,;:]+)$")
SPLIT = re.compile(r"\s*(?:[,;/]|\band\b)\s*")
# same priority as the voucher step (g2t.voucher): specimen_voucher > isolate > culture_collection > clone > strain
VOUCHER_FIELDS = ["specimen_voucher", "isolate", "culture_collection", "clone", "strain"]
MONTHS = {m: i + 1 for i, m in enumerate(["jan", "feb", "mar", "apr", "may", "jun",
                                           "jul", "aug", "sep", "oct", "nov", "dec"])}
RANK = {"high": 2, "medium": 1}


@dataclass
class ReconcileConfig:
    key_column: str = "species_voucher_new"
    organism_column: str = "organism"
    min_confidence: str = "medium"        # "high" = only strong evidence merges
    min_core_length: int = 4


@dataclass
class Evidence:
    strong: list[str] = field(default_factory=list)
    moderate: list[str] = field(default_factory=list)
    conflicts: list[str] = field(default_factory=list)

    def text(self) -> str:
        parts = [f"strong: {'; '.join(self.strong)}" if self.strong else "",
                 f"moderate: {'; '.join(self.moderate)}" if self.moderate else "",
                 f"conflict: {'; '.join(self.conflicts)}" if self.conflicts else ""]
        return " | ".join(p for p in parts if p)


# --------------------------------------------------------------------------- parsing helpers

def _core_of(part: str, min_len: int) -> str:
    v = part
    for _ in range(3):
        v = TRAILING.sub("", v)
        v = GENE_PREFIX.sub("", v)
        v = GENE_SUFFIX.sub("", v)
    m = CORE.search(" " + v)
    if not m:
        return ""
    core = (m.group(1) + re.sub(r"[\s_.\-]", "", m.group(2))).upper()
    return core if len(core) >= min_len else ""


def voucher_cores(raw: str, min_len: int = 4) -> list[str]:
    """Candidate keys of a voucher string: one core per listed identifier (gene names, trailing
    notes such as '(holotype)' and punctuation removed; a year-number pair such as 2014-1234 is kept
    whole), every letters+digits identifier inside a compound voucher (ZMMU_MSU_WS14906_XZ5507 -> XZ5507,
    WS14906), plus a digits-only key ('#123456', >= 5 digits) so that 'USNM 123456' meets '123456'."""
    v = str(raw or "").strip()
    if not v or v.lower() == "nan":
        return []
    out: list[str] = []
    for part in SPLIT.split(v):
        c = _core_of(part.strip(), min_len)
        if c and c not in out:
            out.append(c)
        # compound vouchers carry several identifiers (ZMMU_MSU_WS14906_XZ5507): each one is a candidate key
        for tok in re.split(r"[\s_:\-]+", part.strip()):
            if ID_TOKEN.fullmatch(tok) and not GENE_FULL.fullmatch(tok):
                t = tok.upper()
                if len(t) >= min_len and t not in out:
                    out.append(t)
    for c in list(out):
        digits = re.sub(r"^[A-Z]+", "", c)
        if digits.isdigit() and len(digits) >= 5 and f"#{digits}" not in out:
            out.append(f"#{digits}")
    return out


def voucher_core(raw: str, min_len: int = 4) -> str:
    """Primary core identifier of a voucher string ('' if none)."""
    cores = voucher_cores(raw, min_len)
    return cores[0] if cores else ""


def _genes(df: pd.DataFrame) -> set:
    return {g.strip() for v in _vals(df, "gene_type") for g in v.split(",") if g.strip()}


def _vals(df: pd.DataFrame, col: str) -> set:
    if col not in df:
        return set()
    return {str(x).strip() for x in df[col] if str(x).strip() and str(x).strip().lower() != "nan"}


def _titles(df: pd.DataFrame) -> set:
    out = set()
    for c in df.columns:
        if re.fullmatch(r"Ref\d+Title", str(c)):
            for t in _vals(df, c):
                t = re.sub(r"\W+", " ", t.lower()).strip()
                if t and t != "direct submission":
                    out.add(t)
    return out


def _parse_date(s: str) -> tuple[int | None, int | None, int | None]:
    s = s.strip().lower()
    m = re.fullmatch(r"(\d{4})-(\d{1,2})(?:-(\d{1,2}))?", s)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3)) if m.group(3) else None
    m = re.fullmatch(r"(?:(\d{1,2})-)?([a-z]{3})-(\d{4})", s)
    if m:
        return int(m.group(3)), MONTHS.get(m.group(2)), int(m.group(1)) if m.group(1) else None
    m = re.search(r"(\d{4})", s)
    return (int(m.group(1)), None, None) if m else (None, None, None)


def _dates_compatible(a: str, b: str) -> tuple[bool, bool]:
    """(compatible, identical-to-the-day)."""
    if "/" in a or "/" in b:          # ranges: compare years only
        ya, yb = set(re.findall(r"\d{4}", a)), set(re.findall(r"\d{4}", b))
        return bool(ya & yb) or not (ya and yb), False
    pa, pb = _parse_date(a), _parse_date(b)
    for x, y in zip(pa, pb):
        if x is not None and y is not None and x != y:
            return False, False
    return True, all(x is not None for x in pa) and pa == pb


def _latlon(s: str) -> tuple[float, float] | None:
    m = re.match(r"\s*(-?\d+(?:\.\d+)?)\s*([NS])\s+(-?\d+(?:\.\d+)?)\s*([EW])", s, re.I)
    if not m:
        return None
    lat = float(m.group(1)) * (-1 if m.group(2).upper() == "S" else 1)
    lon = float(m.group(3)) * (-1 if m.group(4).upper() == "W" else 1)
    return lat, lon


def _country(row_vals: set) -> set:
    return {v.split(":")[0].strip().lower() for v in row_vals if v.split(":")[0].strip()}


def _first_author(df: pd.DataFrame) -> set:
    return {a.split(",")[0].strip().lower() for a in _vals(df, "Ref1Authors") if a.split(",")[0].strip()}


# --------------------------------------------------------------------------- evidence

def compare_groups(a: pd.DataFrame, b: pd.DataFrame) -> Evidence:
    ev = Evidence()
    shared = _titles(a) & _titles(b)
    if shared:
        ev.strong.append(f"same publication ('{sorted(shared)[0][:60]}')")
    # collection date
    da, db = _vals(a, "collection_date"), _vals(b, "collection_date")
    if da and db:
        pairs = [_dates_compatible(x, y) for x in da for y in db]
        if not any(c for c, _ in pairs):
            ev.conflicts.append(f"collection dates differ ({', '.join(sorted(da))} vs {', '.join(sorted(db))})")
        elif any(same for _, same in pairs):
            ev.strong.append("same collection date")
    # coordinates
    la = [p for p in map(_latlon, _vals(a, "lat_lon")) if p]
    lb = [p for p in map(_latlon, _vals(b, "lat_lon")) if p]
    if la and lb:
        def dist(x, y):
            dlon = abs(x[1] - y[1]) % 360
            return math.hypot(x[0] - y[0], min(dlon, 360 - dlon))
        d = min(dist(x, y) for x in la for y in lb)
        if d <= 0.01:
            ev.strong.append("same coordinates")
        elif d <= 0.5:
            ev.moderate.append(f"coordinates {d:.2f} deg apart")
        else:
            ev.conflicts.append(f"coordinates {d:.1f} deg apart")
    # country / locality
    ga = _vals(a, "geo_loc_name") | _vals(a, "country")
    gb = _vals(b, "geo_loc_name") | _vals(b, "country")
    ca, cb = _country(ga), _country(gb)
    if ca and cb and not (ca & cb):
        ev.conflicts.append(f"different countries ({', '.join(sorted(ca))} vs {', '.join(sorted(cb))})")
    detailed = {g.lower() for g in ga if ":" in g} & {g.lower() for g in gb if ":" in g}
    if detailed:
        ev.moderate.append("same locality")
    if _vals(a, "collected_by") & _vals(b, "collected_by"):
        ev.moderate.append("same collector")
    if not shared and (_first_author(a) & _first_author(b)):
        ev.moderate.append("same first author")
    return ev


def _decide(ev: Evidence, gene_overlap: bool, min_conf: str) -> tuple[str, str]:
    """-> (confidence or '', decision text)."""
    if ev.conflicts:
        return "", "not merged: conflicting metadata"
    conf = "high" if ev.strong else ("medium" if len(ev.moderate) >= 2 else "")
    if gene_overlap and conf != "high":
        return "", "not merged: same gene in both groups, strong evidence required"
    if not conf:
        return "", "not merged: insufficient evidence"
    if RANK[conf] < RANK[min_conf]:
        return "", f"not merged: {conf} confidence below threshold ({min_conf})"
    return conf, f"merged ({conf})"


# --------------------------------------------------------------------------- main logic

def reconcile_dataframe(df: pd.DataFrame, config: ReconcileConfig | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    config = config or ReconcileConfig()
    df = df.copy()
    key, org = config.key_column, config.organism_column
    if org not in df.columns:
        raise ValueError(f"organism column '{org}' not found: cannot tell organisms apart, refusing to reconcile")
    df[key] = df[key].fillna("").astype(str)
    df["species_voucher_g2t"] = df[key]

    def raw_voucher(row) -> str:
        for c in VOUCHER_FIELDS:
            v = str(row.get(c, "") or "").strip()
            if v and v.lower() != "nan":
                return v
        return ""

    raw = df.apply(raw_voucher, axis=1)
    all_cores = raw.map(lambda v: voucher_cores(v, config.min_core_length))
    df["voucher_core"] = all_cores.map(lambda c: next((x for x in c if not x.startswith("#")), c[0] if c else ""))
    df["match_basis"] = "exact"
    df["match_confidence"] = ""
    df["match_evidence"] = ""

    report: list[dict] = []
    groups = {k: g for k, g in df.groupby(key, sort=False)}
    cores: dict[str, list[str]] = {}
    for k, g in groups.items():
        for c in sorted({c for lst in all_cores[g.index] for c in lst}):
            cores.setdefault(c, []).append(k)

    cache: dict[frozenset, Evidence] = {}

    def evidence(ka: str, kb: str) -> Evidence:
        k = frozenset((ka, kb))
        if k not in cache:
            cache[k] = compare_groups(groups[ka], groups[kb])
        return cache[k]

    def row_for(core, ka, kb):
        ga, gb = groups[ka], groups[kb]
        return dict(voucher_core=core, group_a=ka, group_b=kb,
                    organism_a="; ".join(sorted(_vals(ga, org))), organism_b="; ".join(sorted(_vals(gb, org))),
                    accessions_a="; ".join(map(str, ga["ACCESSION"])) if "ACCESSION" in ga else "",
                    accessions_b="; ".join(map(str, gb["ACCESSION"])) if "ACCESSION" in gb else "",
                    genes_a="; ".join(sorted(_genes(ga))), genes_b="; ".join(sorted(_genes(gb))))

    edges = []
    seen_pairs: set[frozenset] = set()
    for core, keys in cores.items():
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                ka, kb = keys[i], keys[j]
                if frozenset((ka, kb)) in seen_pairs:
                    continue
                seen_pairs.add(frozenset((ka, kb)))
                oa, ob = _vals(groups[ka], org), _vals(groups[kb], org)
                base = row_for(core, ka, kb)
                if not oa or not ob:
                    report.append({**base, "evidence": "", "confidence": "",
                                   "decision": "not merged: organism missing"})
                    continue
                if oa != ob:
                    report.append({**base, "evidence": "", "confidence": "",
                                   "decision": "not merged: different organisms"})
                    continue
                ev = evidence(ka, kb)
                overlap = bool(_genes(groups[ka]) & _genes(groups[kb]))
                conf, decision = _decide(ev, overlap, config.min_confidence)
                report.append({**base, "evidence": ev.text(), "confidence": conf, "decision": decision})
                if conf:
                    edges.append((RANK[conf], len(ev.strong), len(ev.moderate), ka, kb, conf, ev.text()))

    parent = {k: k for k in groups}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    members: dict[str, set[str]] = {k: {k} for k in groups}

    def blocker(ra: str, rb: str) -> str:
        """Check every pair across the two clusters, not only the pairs that shared a core."""
        for x in sorted(members[ra]):
            for y in sorted(members[rb]):
                if _vals(groups[x], org) != _vals(groups[y], org):
                    return "not merged: would join different organisms"
                ev = evidence(x, y)
                if ev.conflicts:
                    return f"not merged: would join conflicting groups ({x} vs {y}: {'; '.join(ev.conflicts)})"
                if _genes(groups[x]) & _genes(groups[y]) and not ev.strong:
                    return (f"not merged: would put the same gene twice in one specimen without strong "
                            f"evidence ({x} vs {y})")
        return ""

    edge_info: dict[str, list[tuple[str, str]]] = {}
    for _, _, _, ka, kb, conf, text in sorted(edges, key=lambda e: (-e[0], -e[1], -e[2], e[3], e[4])):
        ra, rb = find(ka), find(kb)
        if ra == rb:
            continue
        why = blocker(ra, rb)
        if why:
            for r in report:
                if {r["group_a"], r["group_b"]} == {ka, kb}:
                    r["decision"] = why
            continue
        parent[rb] = ra
        members[ra] |= members.pop(rb)
        edge_info.setdefault(ra, []).extend(edge_info.pop(rb, []))
        edge_info[ra].append((conf, f"{ka} + {kb}: {text}"))

    used_keys = set(groups)
    for root in sorted(members):
        member_keys = members[root]
        if len(member_keys) < 2:
            continue
        sub = df[df[key].isin(member_keys)]
        orgs = sorted(_vals(sub, org))
        organism = orgs[0] if orgs else "unknown"
        core = sorted(set(sub["voucher_core"]) - {""})[0]
        base_key = f"{re.sub(r'[^A-Za-z0-9]+', '_', organism).strip('_')}_{core.lstrip('#')}"
        new_key, n = base_key, 1
        while new_key in used_keys - member_keys:
            new_key, n = f"{base_key}_reconciled{n}", n + 1
        used_keys.add(new_key)
        confs = [c for c, _ in edge_info.get(root, [])]
        conf = "high" if confs and all(c == "high" for c in confs) else "medium"
        mask = df[key].isin(member_keys)
        df.loc[mask, key] = new_key
        df.loc[mask, "match_basis"] = "reconciled"
        df.loc[mask, "match_confidence"] = conf
        df.loc[mask, "match_evidence"] = " || ".join(t for _, t in edge_info.get(root, []))

    cols = ["voucher_core", "group_a", "group_b", "organism_a", "organism_b", "genes_a", "genes_b",
            "accessions_a", "accessions_b", "evidence", "confidence", "decision"]
    return df, pd.DataFrame(report, columns=cols)


def reconcile(input_file: str, output_dir: str, config: ReconcileConfig | None = None,
              output_name: str = "reconciled_species_voucher.csv") -> StepResult:
    """Run Step 3b on the voucher-step output."""
    t0 = time.time()
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = read_table(input_file)
    out, report = reconcile_dataframe(df, config)
    out_file = out_dir / output_name
    out.to_csv(out_file, index=False)
    report.to_csv(out_dir / "reconcile_report.csv", index=False)
    return StepResult(success=True, output_file=str(out_file), rows=len(out), elapsed=time.time() - t0)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Evidence-based merging of voucher variants (Step 3b)")
    p.add_argument("-i", "--input", required=True, help="updated_species_voucher.csv from g2t-voucher")
    p.add_argument("-o", "--output_dir", required=True)
    p.add_argument("--min_confidence", choices=["high", "medium"], default="medium",
                   help="high = merge only with strong evidence (shared paper/date/coordinates)")
    a = p.parse_args(argv)
    r = reconcile(a.input, a.output_dir, ReconcileConfig(min_confidence=a.min_confidence))
    rep = pd.read_csv(Path(a.output_dir) / "reconcile_report.csv")
    merged = rep["decision"].str.startswith("merged").sum() if len(rep) else 0
    print(f"Reconciled: {r.output_file} ({r.rows} rows); candidate pairs {len(rep)}, merged {merged}")


if __name__ == "__main__":
    main()
