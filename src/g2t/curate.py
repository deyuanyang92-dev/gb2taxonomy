#!/usr/bin/env python3
"""
Voucher standardisation and metadata curation (Step 5)
======================================================

1. ``standardize_vouchers``: one voucher per specimen. GenBank submitters write the same voucher
   differently per record (``COI_ZMMU_MSU_WS2585`` / ``28S_ZMMU_WS2585`` / ``WS2585``,
   ``ZMMU:WS30980`` / ``ZMMU_WS30980``). Each written form is normalised (gene names removed,
   separators unified), the most complete form is kept and written in the INSDC format
   ``[<institution-code>:[<collection-code>:]]<specimen_id>`` (``ZMMU:MSU:WS2585``, ``ZMMU:WS30980``).
   The organize step writes it to ``voucher_standardized`` next to ``voucher_as_submitted``
   (all forms as found in GenBank).

2. ``apply_updates`` / ``g2t-curate``: GenBank metadata is often outdated (coordinates, species names,
   publications added after submission). A user table with one row per correction updates the
   matrix. Rows are located by ``accession`` (any gene column; version optional) and/or ``voucher``
   (any written form; optionally narrowed with ``organism_match``). Empty cells change nothing.
   Every change is logged (old -> new, matched by); rows that match nothing, several rows, or
   contradict each other are reported and not applied. The input matrix is never modified.

CLI:
    g2t-curate -m organized_species_voucher.csv --template corrections.xlsx   # keeps a baseline copy
    g2t-curate -m organized_species_voucher.csv -u corrections.csv -o curated_matrix.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

GENE_TOKENS = (r"COI|CO1|COX1|COII|COX2|COIII|COX3|CYTB|COB|ND\d|NAD\d|ATP\d|12S|16S|18S|28S|5\.8S|"
               r"H3|ITS1?|ITS2|EF1A?|EF1ALPHA|RRNL|RRNS|LSU|SSU|MTDNA|MITO")
_PREFIX = re.compile(rf"^(?:{GENE_TOKENS})[\s_\-:]+", re.I)
_SUFFIX = re.compile(rf"[\s_\-:]+(?:{GENE_TOKENS})$", re.I)
KEY_COLUMNS = ["row_id", "voucher", "accession", "organism_match"]
BASELINE_SHEET = "baseline (do not edit)"
TEMPLATE_FIELDS = ["organism", "lat_lon", "geo_loc_name", "country", "collection_date", "collected_by",
                   "identified_by", "Ref1Authors", "Ref1Title", "Ref1Journal", "curation_source"]
DEFAULT_GENES = ["mtgenome", "coi", "16s", "12s", "cob", "cox2", "cox3", "18s", "28s", "its1-its2",
                 "18-28s", "ef-1", "h3"]


# --------------------------------------------------------------------------- vouchers

def normalize_voucher(raw) -> str:
    """'COI_ZMMU_MSU_WS2585' -> 'ZMMU_MSU_WS2585'; 'ZMMU:WS30980' -> 'ZMMU_WS30980'."""
    v = "" if raw is None or (isinstance(raw, float) and pd.isna(raw)) else str(raw).strip()
    if not v or v.lower() == "nan":
        return ""
    for _ in range(2):
        v = _PREFIX.sub("", v)
        v = _SUFFIX.sub("", v)
    v = re.sub(r"[\s:]+", "_", v)
    v = re.sub(r"_+", "_", v).strip("_")
    return v


def _tokens(v: str) -> list[str]:
    return [t for t in re.split(r"[_\-]", v.upper()) if t]


def _contains(big: str, small: str) -> bool:
    """small's tokens appear in big, in order (e.g. WS2585 in ZMMU_MSU_WS2585)."""
    b, s = _tokens(big), _tokens(small)
    i = 0
    for t in b:
        if i < len(s) and t == s[i]:
            i += 1
    return i == len(s)


def to_insdc(norm: str) -> str:
    """Normalised voucher -> INSDC /specimen_voucher form "[<institution-code>:[<collection-code>:]]<specimen_id>".
    Leading letters-only tokens are the codes, the rest is the specimen id:
    ZMMU_WS30980 -> ZMMU:WS30980; ZMMU_MSU_WS12387_XZ5022 -> ZMMU:MSU:WS12387_XZ5022; WS0397 -> WS0397.
    With more than two code tokens the form is left unchanged (cannot tell institution from collection)."""
    if not norm:
        return ""
    toks = norm.split("_")
    codes = []
    for t in toks[:-1]:
        if re.fullmatch(r"[A-Za-z]{2,10}", t):
            codes.append(t)
        else:
            break
    if not codes or len(codes) > 2:
        return norm
    return ":".join(codes + ["_".join(toks[len(codes):])])


def standardize_vouchers(values) -> tuple[str, bool]:
    """-> (standardised voucher in INSDC form, consistent). The most complete normalised form is kept
    (see to_insdc); consistent is False when some form is not contained in it (more than a
    prefix/separator difference)."""
    forms: list[str] = []
    for raw in values:
        for part in str(raw if raw is not None else "").split(";"):
            n = normalize_voucher(part)
            if n and n not in forms:
                forms.append(n)
    if not forms:
        return "", True
    best = max(forms, key=lambda f: (len(_tokens(f)), len(f), -forms.index(f)))
    return to_insdc(best), all(_contains(best, f) for f in forms)


# --------------------------------------------------------------------------- curation

def _acc(a) -> str:
    return str(a).strip().split(".")[0].upper()


def _gene_cols(matrix: pd.DataFrame, gene_columns=None) -> list[str]:
    return [g for g in (gene_columns or DEFAULT_GENES) if g in matrix.columns]


def _voucher_forms(row: pd.Series) -> set[str]:
    forms = set()
    for c in ("voucher_standardized", "voucher_as_submitted", "specimen_voucher", "isolate"):
        if c in row and str(row[c]).strip() and str(row[c]).lower() != "nan":
            for part in str(row[c]).split(";"):
                n = normalize_voucher(part).upper()
                if n:
                    forms.add(n)
    return forms


def _blank(v) -> bool:
    return v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() == "" or str(v).lower() == "nan"


def apply_updates(matrix: pd.DataFrame, updates: pd.DataFrame, gene_columns=None,
                  key_column: str = "species_voucher_new",
                  baseline: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """-> (curated matrix, change log, problems). See module docstring.

    ``baseline``: the values a pre-filled template was made from (rows joined on ``row_id``). With it,
    only cells the user edited count as corrections, so unedited template cells never overwrite values
    that GenBank updated later; an edited cell whose GenBank value also changed is reported, not applied."""
    base_map: dict[str, dict] = {}
    if baseline is not None and "row_id" in baseline.columns:
        base_map = {str(r["row_id"]).strip(): r for _, r in baseline.fillna("").astype(str).iterrows()}
    out = matrix.copy().astype(object)
    genes = _gene_cols(out, gene_columns)
    acc_index: dict[str, set[int]] = {}
    for i, row in out.iterrows():
        for g in genes:
            for a in str(row[g]).split(";"):
                if a.strip() and a.strip().lower() != "nan":
                    acc_index.setdefault(_acc(a), set()).add(i)
    forms = {i: _voucher_forms(row) for i, row in out.iterrows()}
    if "curated_fields" not in out.columns:
        out["curated_fields"] = ""
    log, problems = [], []
    upd = updates.copy()
    upd.columns = [str(c).strip() for c in upd.columns]
    fields = [c for c in upd.columns if c not in KEY_COLUMNS]
    for c in fields:
        if c not in out.columns:
            out[c] = ""
            if c not in TEMPLATE_FIELDS:
                problems.append({"update_row": "", "keys": "", "problem": f"new column '{c}' added to the matrix"})
    for n, r in upd.iterrows():
        acc = "" if _blank(r.get("accession")) else str(r.get("accession")).strip()
        vou = "" if _blank(r.get("voucher")) else str(r.get("voucher")).strip()
        org = "" if _blank(r.get("organism_match")) else str(r.get("organism_match")).strip()
        keys = "; ".join(x for x in (f"accession={acc}" if acc else "", f"voucher={vou}" if vou else "",
                                     f"organism_match={org}" if org else "") if x)
        rownum = n + 2                                   # spreadsheet row (header = 1)
        if not acc and not vou:
            problems.append({"update_row": rownum, "keys": keys, "problem": "no accession or voucher given"})
            continue
        hit_a = acc_index.get(_acc(acc), set()) if acc else None
        q = normalize_voucher(vou).upper()
        hit_v = {i for i, f in forms.items() if any(_contains(x, q) or _contains(q, x) for x in f)} if vou else None
        if org:
            keep = {i for i in out.index if str(out.at[i, "organism"]).strip() == org} if "organism" in out else set()
            hit_a = hit_a & keep if hit_a is not None else None
            hit_v = hit_v & keep if hit_v is not None else None
        if hit_a is not None and hit_v is not None and hit_a != hit_v:
            problems.append({"update_row": rownum, "keys": keys,
                             "problem": "accession and voucher point to different rows; not applied"})
            continue
        hits = hit_a if hit_a is not None else hit_v
        if not hits:
            problems.append({"update_row": rownum, "keys": keys, "problem": "no matrix row matches; not applied"})
            continue
        if len(hits) > 1:
            names = ", ".join(str(out.at[i, key_column]) for i in sorted(hits)) if key_column in out else ""
            problems.append({"update_row": rownum, "keys": keys,
                             "problem": f"matches {len(hits)} matrix rows ({names}); add organism_match or "
                                        "use an accession; not applied"})
            continue
        i = next(iter(hits))
        by = "accession" if acc else "voucher"
        base = base_map.get(str(r.get("row_id", "")).strip()) if base_map else None
        for c in fields:
            new = r[c]
            if _blank(new):
                continue
            old = out.at[i, c]
            new = str(new).strip()
            cur = "" if _blank(old) else str(old).strip()
            if base is not None and c in base.index:
                was = str(base[c]).strip()
                if new == was:
                    continue                      # not edited by the user
                if cur != was:
                    problems.append({"update_row": rownum, "keys": keys,
                                     "problem": f"{c}: GenBank value changed since the template was made "
                                                f"('{was}' -> '{cur}'); your value '{new}' not applied"})
                    continue
            if cur == new:
                continue
            out.at[i, c] = new
            done = [x for x in str(out.at[i, "curated_fields"]).split("; ") if x]
            if c not in done:
                out.at[i, "curated_fields"] = "; ".join(done + [c])
            log.append({"specimen": out.at[i, key_column] if key_column in out else i, "field": c,
                        "old_value": "" if _blank(old) else old, "new_value": new, "matched_by": by,
                        "keys": keys, "update_row": rownum})
    log_df = pd.DataFrame(log, columns=["specimen", "field", "old_value", "new_value", "matched_by", "keys",
                                        "update_row"])
    prob_df = pd.DataFrame(problems, columns=["update_row", "keys", "problem"])
    return out, log_df, prob_df


def make_template(matrix: pd.DataFrame, gene_columns=None, fields=None) -> pd.DataFrame:
    """One row per specimen, pre-filled with the current values, for the user to correct."""
    genes = _gene_cols(matrix, gene_columns)
    rows = []
    for _, r in matrix.iterrows():
        acc = next((str(r[g]).split(";")[0].strip() for g in genes
                    if str(r[g]).strip() and str(r[g]).lower() != "nan"), "")
        voucher = r.get("voucher_standardized", "")
        row = {"row_id": len(rows) + 1, "voucher": "" if _blank(voucher) else voucher, "accession": acc,
               "organism_match": r.get("organism", "")}
        for f in fields or TEMPLATE_FIELDS:
            v = r.get(f, "")
            row[f] = "" if _blank(v) else v
        rows.append(row)
    return pd.DataFrame(rows)


def _read(path: str, sheet=0) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p, dtype=str, sheet_name=sheet).fillna("")
    sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    return pd.read_csv(p, dtype=str, sep=sep, keep_default_na=False)


def write_template(template: pd.DataFrame, path: str) -> None:
    """xlsx: sheets 'corrections' + baseline; csv: <name>.csv + <name>.baseline.csv."""
    p = Path(path)
    if p.suffix.lower() == ".xlsx":
        with pd.ExcelWriter(p) as w:
            template.to_excel(w, sheet_name="corrections", index=False)
            template.to_excel(w, sheet_name=BASELINE_SHEET, index=False)
    else:
        template.to_csv(p, index=False)
        template.to_csv(p.with_name(p.stem + ".baseline.csv"), index=False)


def read_corrections(path: str) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """-> (corrections, baseline or None)."""
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        sheets = pd.read_excel(p, dtype=str, sheet_name=None)
        names = list(sheets)
        upd = sheets["corrections"] if "corrections" in sheets else sheets[names[0]]
        base = sheets.get(BASELINE_SHEET)
        return upd.fillna(""), (base.fillna("") if base is not None else None)
    b = p.with_name(p.stem + ".baseline.csv")
    return _read(path), (_read(str(b)) if b.exists() else None)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Correct matrix metadata from a table keyed by accession or voucher")
    p.add_argument("-m", "--matrix", required=True, help="organized_species_voucher.csv (or .xlsx)")
    p.add_argument("-u", "--updates", help="corrections table (.csv/.tsv/.xlsx): accession and/or voucher "
                                           "[+ organism_match] + columns to set")
    p.add_argument("-o", "--output", help="curated matrix (.csv); log and problems are written next to it")
    p.add_argument("--template", help="write a pre-filled corrections template (.csv or .xlsx) and exit")
    a = p.parse_args(argv)
    matrix = _read(a.matrix)
    if a.template:
        t = make_template(matrix)
        write_template(t, a.template)
        print(f"Template: {a.template} ({len(t)} specimens). Edit the cells to correct, leave blank to keep.")
        return
    if not (a.updates and a.output):
        p.error("-u/--updates and -o/--output are required (or use --template)")
    upd, base = read_corrections(a.updates)
    out, log, problems = apply_updates(matrix, upd, baseline=base)
    o = Path(a.output)
    out.to_csv(o, index=False)
    log.to_csv(o.with_name(o.stem + "_curation_log.csv"), index=False)
    problems.to_csv(o.with_name(o.stem + "_curation_problems.csv"), index=False)
    print(f"Curated: {o} ({len(log)} changes in {log['specimen'].nunique() if len(log) else 0} specimens; "
          f"{len(problems)} problems)")


if __name__ == "__main__":
    main()
