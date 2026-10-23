#!/usr/bin/env python3
"""
Voucher Reconciler (Step 3b)
============================
Merges specimen groups whose vouchers are written differently but refer to the
same specimen, e.g. ``COI_ZMMU_MSU_WS399`` / ``28S_ZMMU_WS399`` / ``WS399``.

Exact string matching (Step 3) splits such specimens into one row per gene.
This step decides merges from evidence in the GenBank records themselves,
never from the voucher string alone:

1. Candidates: groups whose identifiers share a core. Identifiers are taken
   from every field (specimen_voucher, isolate, culture_collection, clone,
   strain); a core is what remains after gene names, collection acronyms and
   separators are removed (>= 4 characters; COI_ZMMU_MSU_WS399 -> WS399).
   Pairs must be the same organism, or the same genus with an identical
   identifier (ignoring case and punctuation): one specimen deposited under two
   spellings of its name (Bathykurila sp. / Bathykurila sp. AH-2023, both
   SIO:BIC:A10920). Other pairs are reported, never merged; an identical
   identifier in different genera is also flagged in match_evidence.
2. Evidence for each candidate pair:
     strong    - identical identifier (ignoring case and punctuation);
                 a shared publication title (not "Direct Submission");
                 identical full collection date; coordinates within 0.01 deg
     moderate  - same first author; same collector; same detailed locality;
                 coordinates within 0.5 deg; an identical identifier when both
                 groups already contain the same gene
     conflict  - the shared core is followed by different endings (lot or
                 sub-sample numbers: BIOUG06836-E01 / -E02, SIO:BIC:A9392-1 /
                 -2, CE438 / CE438_2); incompatible collection dates; different
                 countries; coordinates > 0.5 deg apart
3. Decision: any conflict -> not merged; >= 1 strong -> merge ("high");
   >= 2 moderate -> merge ("medium"); otherwise not merged. If both groups
   already contain the same gene, only the same collection date or the same
   coordinates can merge them (a paper reports many specimens); such pairs left
   apart are noted in match_evidence ("similar voucher ... check").
   Merges are applied strongest first, and a merge is skipped if any two groups
   it would join conflict, differ in organism (other than the genus rule
   above) or share a gene without date/coordinate evidence.

Outputs (in ``output_dir``):
  reconciled_species_voucher.csv  input rows + species_voucher_g2t, voucher_core,
                                  match_basis (exact/reconciled),
                                  match_confidence, match_evidence;
                                  species_voucher_new holds the reconciled key
  reconcile_report.csv            one row per candidate pair with the decision
"""

from __future__ import annotations

import argparse
import logging
import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from g2t.utils import StepResult, read_table

logger = logging.getLogger(__name__)

GENE_TOKENS = (
    r"COI|CO1|COX1|COII|COX2|COIII|COX3|CYTB|COB|ND\d|NAD\d|ATP\d|12S|16S|18S|28S|5\.8S|"
    r"H3|ITS1?|ITS2|EF1A?|EF1ALPHA|RRNL|RRNS|LSU|SSU|MTDNA|MITO"
)
GENE_PREFIX = re.compile(rf"^(?:{GENE_TOKENS})[\s_\-:]+", re.I)
GENE_SUFFIX = re.compile(rf"[\s_\-:]+(?:{GENE_TOKENS})$", re.I)
CORE = re.compile(r"(?:^|[^A-Za-z0-9])([A-Za-z]{0,6})[\s_\-]?((?:(?:1[89]|20)\d{2}[\s_.\-])?\d{2,8}[A-Za-z]?)$")
ID_TOKEN = re.compile(r"[A-Za-z]{1,6}\d{3,8}[A-Za-z]?")
GENE_FULL = re.compile(rf"(?:{GENE_TOKENS})", re.I)
TRAILING = re.compile(r"(\s*\([^()]*\)|[\s.,;:]+)$")
SPLIT = re.compile(r"\s*(?:[,;/]|\band\b)\s*")
# same priority as the voucher step (g2t.voucher): specimen_voucher > isolate > culture_collection > clone > strain
VOUCHER_FIELDS = ["specimen_voucher", "isolate", "culture_collection", "clone", "strain"]
MONTHS = {
    m: i + 1 for i, m in enumerate(["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"])
}
RANK = {"high": 2, "medium": 1}


@dataclass
class ReconcileConfig:
    key_column: str = "species_voucher_new"
    organism_column: str = "organism"
    min_confidence: str = "medium"  # "high" = only strong evidence merges
    min_core_length: int = 4


@dataclass
class Evidence:
    strong: list[str] = field(default_factory=list)
    moderate: list[str] = field(default_factory=list)
    conflicts: list[str] = field(default_factory=list)

    def text(self) -> str:
        parts = [
            f"strong: {'; '.join(self.strong)}" if self.strong else "",
            f"moderate: {'; '.join(self.moderate)}" if self.moderate else "",
            f"conflict: {'; '.join(self.conflicts)}" if self.conflicts else "",
        ]
        return " | ".join(p for p in parts if p)


# --------------------------------------------------------------------------- parsing helpers


def _clean_id(part: str) -> str:
    """Identifier without gene-name prefix/suffix and trailing notes ('COI_WS399 (holotype)' -> 'WS399')."""
    v = part
    for _ in range(3):
        v = TRAILING.sub("", v)
        v = GENE_PREFIX.sub("", v)
        v = GENE_SUFFIX.sub("", v)
    return v


def _core_of(part: str, min_len: int) -> str:
    v = _clean_id(part)
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


def _norm_voucher(v: str) -> str:
    """'SIO:BIC:A10920', 'SIO-BIC A10920' -> 'SIOBICA10920'."""
    return re.sub(r"[^A-Z0-9]", "", str(v).upper())


def _sub_ids(ids: list[str], core: str) -> set[str]:
    """What follows the core in each identifier, after normalisation: 'SIO:BIC:A9392-1' + core A9392 -> '1';
    'COI_ZMMU_WS399' + core WS399 -> ''. Identifiers that do not contain the core are skipped."""
    c = core.lstrip("#")
    out = set()
    for v in ids:
        n = _norm_voucher(_clean_id(v))
        i = n.rfind(c)
        if i >= 0:
            end = n[i + len(c) :]
            # a second full identifier is not a sub-number (ZMMU_MSU_WS14906_XZ5507 also names WS14906)
            out.add("" if ID_TOKEN.fullmatch(end) and len(end) >= 4 else end)
    return out


def _genus(organism: str) -> str:
    """First word of the name after open-nomenclature qualifiers ('cf. Alluroididae sp.' -> 'alluroididae')."""
    words = [w for w in str(organism).split() if w.lower().rstrip(".") not in {"cf", "aff", "nr", "?"}]
    return words[0].lower() if words else ""


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
    if "/" in a or "/" in b:  # ranges: compare years only
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


def _specimen_level(strong: list[str]) -> list[str]:
    """Strong evidence that points to one specimen: a shared paper does not (a paper has many specimens)."""
    return [e for e in strong if not e.startswith("same publication")]


def _decide(ev: Evidence, gene_overlap: bool, min_conf: str) -> tuple[str, str]:
    """-> (confidence or '', decision text)."""
    if ev.conflicts:
        return "", f"not merged: conflicting metadata ({'; '.join(ev.conflicts)})"
    conf = "high" if ev.strong else ("medium" if len(ev.moderate) >= 2 else "")
    if gene_overlap and not _specimen_level(ev.strong):
        return "", "not merged: same gene in both groups, same date/coordinates required"
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
            if any(ch.isalnum() for ch in v) and v.lower() != "nan":  # same rule as voucher.py (L-22)
                return v
        return ""

    raw = df.apply(raw_voucher, axis=1)
    raw_cores = raw.map(lambda v: voucher_cores(v, config.min_core_length))
    df["voucher_core"] = raw_cores.map(lambda c: next((x for x in c if not x.startswith("#")), c[0] if c else ""))

    # candidates come from every identifier field: one record may give the code as specimen_voucher,
    # another as isolate (GU901728: voucher SMNH<SWE>:108431 + isolate CE438; PP994571: isolate CE438)
    def row_ids(row) -> list[str]:
        out = []
        for c in VOUCHER_FIELDS:
            v = str(row.get(c, "") or "").strip()
            if any(ch.isalnum() for ch in v) and v.lower() != "nan":
                out += [p for p in SPLIT.split(v) if p.strip()]
        return out

    ids = df.apply(row_ids, axis=1)
    all_cores = ids.map(
        lambda lst: list(dict.fromkeys(c for v in lst for c in voucher_cores(v, config.min_core_length)))
    )
    df["match_basis"] = "exact"
    df["match_confidence"] = ""
    df["match_evidence"] = ""

    report: list[dict] = []
    groups = {k: g for k, g in df.groupby(key, sort=False)}
    group_ids = {k: [v for lst in ids[g.index] for v in lst] for k, g in groups.items()}
    group_cores = {k: {c for lst in all_cores[g.index] for c in lst} for k, g in groups.items()}
    # every identifier field, not only the one that made the key: one record may carry the code as
    # specimen_voucher, another as isolate (MF001136 'personal collection:...:CE14346' + isolate CE14346)
    norm_vouchers = {
        k: {
            n
            for c in VOUCHER_FIELDS
            if c in g.columns
            for n in (_norm_voucher(v) for v in g[c].fillna("").astype(str))
            if len(n) >= config.min_core_length and n != "NAN"
        }
        for k, g in groups.items()
    }

    def genera(k: str) -> set:
        return {_genus(o) for o in _vals(groups[k], org)}

    def same_genus(ka: str, kb: str) -> bool:
        return len(genera(ka) | genera(kb)) == 1

    cores: dict[str, list[str]] = {}
    for k in groups:
        for c in sorted(group_cores[k]):
            cores.setdefault(c, []).append(k)

    cache: dict[frozenset, Evidence] = {}

    def evidence(ka: str, kb: str) -> Evidence:
        k = frozenset((ka, kb))
        if k not in cache:
            ev = compare_groups(groups[ka], groups[kb])
            # same core, different endings: lot or sub-sample numbers of different specimens
            # (BIOUG06836-E01 / -E02, SIO:BIC:A9392-1 / -2, CE438 / CE438_2)
            shared = sorted(group_cores[ka] & group_cores[kb])
            if shared:
                ends = [(_sub_ids(group_ids[ka], c), _sub_ids(group_ids[kb], c)) for c in shared]
                ends = [(a, b) for a, b in ends if a and b]
                if ends and all(not (a & b) for a, b in ends):
                    a, b = ends[0]
                    ev.conflicts.append(
                        f"different sub-identifiers after {shared[0].lstrip('#')} "
                        f"({', '.join(sorted(x or '-' for x in a))} vs {', '.join(sorted(x or '-' for x in b))})"
                    )
            if norm_vouchers[ka] & norm_vouchers[kb]:
                oa, ob = _vals(groups[ka], org), _vals(groups[kb], org)
                names = f" (organism names differ: {'; '.join(sorted(oa | ob))})" if oa != ob else ""
                # both groups have the same gene: an identical label alone does not make them one specimen
                if _genes(groups[ka]) & _genes(groups[kb]):
                    ev.moderate.insert(0, "identical voucher" + names)
                else:
                    ev.strong.insert(0, "identical voucher" + names)
            cache[k] = ev
        return cache[k]

    def row_for(core, ka, kb):
        ga, gb = groups[ka], groups[kb]
        return dict(
            voucher_core=core,
            group_a=ka,
            group_b=kb,
            organism_a="; ".join(sorted(_vals(ga, org))),
            organism_b="; ".join(sorted(_vals(gb, org))),
            accessions_a="; ".join(map(str, ga["ACCESSION"])) if "ACCESSION" in ga else "",
            accessions_b="; ".join(map(str, gb["ACCESSION"])) if "ACCESSION" in gb else "",
            genes_a="; ".join(sorted(_genes(ga))),
            genes_b="; ".join(sorted(_genes(gb))),
        )

    edges = []
    flagged: dict[str, set[str]] = {k: set() for k in groups}
    similar: dict[str, set[str]] = {k: set() for k in groups}
    seen_pairs: set[frozenset] = set()
    for core, keys in cores.items():
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                ka, kb = sorted((keys[i], keys[j]))  # edge order independent of input row order
                if frozenset((ka, kb)) in seen_pairs:
                    continue
                seen_pairs.add(frozenset((ka, kb)))
                oa, ob = _vals(groups[ka], org), _vals(groups[kb], org)
                base = row_for(core, ka, kb)
                if not oa or not ob:
                    report.append(
                        {**base, "evidence": "", "confidence": "", "decision": "not merged: organism missing"}
                    )
                    continue
                identical = bool(norm_vouchers[ka] & norm_vouchers[kb])
                if oa != ob and not (identical and same_genus(ka, kb)):
                    if identical:
                        decision = "not merged: same voucher, different genera (check)"
                        flagged[ka].add(kb)
                        flagged[kb].add(ka)
                    else:
                        decision = "not merged: different organisms"
                    report.append({**base, "evidence": "", "confidence": "", "decision": decision})
                    continue
                ev = evidence(ka, kb)
                overlap = bool(_genes(groups[ka]) & _genes(groups[kb]))
                conf, decision = _decide(ev, overlap, config.min_confidence)
                if overlap and not conf and not ev.conflicts and (ev.strong or ev.moderate):
                    similar[ka].add(kb)
                    similar[kb].add(ka)
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
                if _vals(groups[x], org) != _vals(groups[y], org) and not (
                    same_genus(x, y) and norm_vouchers[x] & norm_vouchers[y]
                ):
                    return "not merged: would join different organisms"
                ev = evidence(x, y)
                if ev.conflicts:
                    return f"not merged: would join conflicting groups ({x} vs {y}: {'; '.join(ev.conflicts)})"
                if _genes(groups[x]) & _genes(groups[y]) and not _specimen_level(ev.strong):
                    return (
                        f"not merged: would put the same gene twice in one specimen without strong "
                        f"evidence ({x} vs {y})"
                    )
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
        # the key's core may come from another identifier field than voucher_core ('s.n.' + isolate AB1234)
        shared = set.intersection(*(group_cores[k] for k in member_keys))
        key_cores = (
            sorted(shared - {""})
            or sorted(set(sub["voucher_core"]) - {""})
            or sorted(set().union(*all_cores[sub.index]))
        )
        core = next((c for c in key_cores if not c.startswith("#")), key_cores[0])
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

    def add_note(k: str, note: str) -> None:
        mask = df[key] == k
        df.loc[mask, "match_evidence"] = [f"{e} || {note}" if e else note for e in df.loc[mask, "match_evidence"]]

    n_flag = 0
    for k, others in flagged.items():
        if others:
            n_flag += 1
            add_note(k, f"same voucher as {'; '.join(sorted(others))} (different genus, not merged; check)")
    for k, others in similar.items():
        if others:
            add_note(k, f"similar voucher to {'; '.join(sorted(others))} (both have the same gene, not merged; check)")
    if n_flag:
        logger.warning(
            f"Reconcile: {n_flag} specimen groups share a voucher with a group of another genus and were kept "
            "apart; see match_evidence and reconcile_report.csv"
        )

    cols = [
        "voucher_core",
        "group_a",
        "group_b",
        "organism_a",
        "organism_b",
        "genes_a",
        "genes_b",
        "accessions_a",
        "accessions_b",
        "evidence",
        "confidence",
        "decision",
    ]
    return df, pd.DataFrame(report, columns=cols)


def reconcile(
    input_file: str,
    output_dir: str,
    config: ReconcileConfig | None = None,
    output_name: str = "reconciled_species_voucher.csv",
) -> StepResult:
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
    p.add_argument(
        "--min_confidence",
        choices=["high", "medium"],
        default="medium",
        help="high = merge only with strong evidence (shared paper/date/coordinates)",
    )
    a = p.parse_args(argv)
    r = reconcile(a.input, a.output_dir, ReconcileConfig(min_confidence=a.min_confidence))
    rep = pd.read_csv(Path(a.output_dir) / "reconcile_report.csv")
    merged = rep["decision"].str.startswith("merged").sum() if len(rep) else 0
    print(f"Reconciled: {r.output_file} ({r.rows} rows); candidate pairs {len(rep)}, merged {merged}")


if __name__ == "__main__":
    main()
