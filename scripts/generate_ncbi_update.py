#!/usr/bin/env python3
"""
Generate NCBI Source Modifiers Table for updating GenBank records.

Reads literature metadata and GB current state, produces a tab-delimited
file in NCBI's required format for submission to gb-admin@ncbi.nlm.nih.gov.

Only fills fields that are missing from existing GB records.

Usage:
    python scripts/generate_ncbi_update.py
"""

import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

BASE = Path("/media/deyuan/217ce44c-b5de-45ed-8720-deebdff85ece/CLAUDE/ceshi/ceshi5")
EXCEL = BASE / "zlae141_suppl_supplementary_tables_s1.xlsx"
FINAL_CSV = BASE / "haplotaxis_output3/gb_metadata/final.csv"
OUTPUT = BASE / "ncbi_source_update.tsv"

MONTH_ABBR = {
    "jan": "Jan", "feb": "Feb", "mar": "Mar", "apr": "Apr",
    "may": "May", "jun": "Jun", "jul": "Jul", "aug": "Aug",
    "sep": "Sep", "oct": "Oct", "nov": "Nov", "dec": "Dec",
}


def is_empty(val) -> bool:
    s = str(val).strip()
    return s == "" or s == "nan" or pd.isna(val)


def normalize_date(raw: str) -> str:
    """Convert literature date to NCBI DD-Mon-YYYY format."""
    raw = str(raw).strip()
    if not raw or raw == "nan":
        return ""

    # Handle year-only: "1993"
    if re.match(r"^\d{4}$", raw):
        return raw

    # Handle "Mon YYYY" like "Dec 2000"
    m = re.match(r"^(\w{3})\s+(\d{4})$", raw)
    if m:
        mon = MONTH_ABBR.get(m.group(1).lower(), m.group(1))
        return f"{mon}-{m.group(2)}"

    # Handle "DD-Mon-YYYY" already
    m = re.match(r"^(\d{1,2})-(\w{3})-(\d{4})$", raw)
    if m:
        return raw

    # Handle "DD Mon YYYY" or "D Mon YYYY" (e.g. "08 Nov 2005", "8 Apr 2009")
    m = re.match(r"^(\d{1,2})\s+(\w+\.?)\s+(\d{4})$", raw)
    if m:
        day = m.group(1).zfill(2)
        mon_raw = m.group(2).rstrip(".").lower()
        mon = MONTH_ABBR.get(mon_raw, mon_raw.capitalize())
        return f"{day}-{mon}-{m.group(3)}"

    # Handle date ranges: "18-20 Apr 2012" → use first date
    m = re.match(r"^(\d{1,2})-(\d{1,2})\s+(\w+\.?)\s+(\d{4})$", raw)
    if m:
        day = m.group(1).zfill(2)
        mon_raw = m.group(3).rstrip(".").lower()
        mon = MONTH_ABBR.get(mon_raw, mon_raw.capitalize())
        return f"{day}-{mon}-{m.group(4)}"

    # Handle "Mon-Mon YYYY" like "Jul-Aug 2012"
    m = re.match(r"^(\w{3})-(\w{3})\s+(\d{4})$", raw)
    if m:
        mon = MONTH_ABBR.get(m.group(1).lower(), m.group(1))
        return f"{mon}-{m.group(3)}"

    return raw


def format_geo_loc_name(country, state, municipality, locality) -> str:
    """Build NCBI geo_loc_name: 'Country: State, Municipality, Locality'"""
    parts = []
    if state and str(state).strip() and str(state).strip() != "nan":
        parts.append(str(state).strip())
    if municipality and str(municipality).strip() and str(municipality).strip() != "nan":
        parts.append(str(municipality).strip())
    if locality and str(locality).strip() and str(locality).strip() != "nan":
        parts.append(str(locality).strip())

    country = str(country).strip()
    if parts:
        return f"{country}: {', '.join(parts)}"
    return country


def format_lat_lon(n, e) -> str:
    """Format coordinates as '57.68 N 11.96 E'."""
    try:
        lat = float(str(n).strip())
        lon = float(str(e).strip())
    except (ValueError, TypeError):
        return ""
    lat_dir = "N" if lat >= 0 else "S"
    lon_dir = "E" if lon >= 0 else "W"
    return f"{abs(lat):.4f} {lat_dir} {abs(lon):.4f} {lon_dir}"


def format_voucher(voucher: str) -> str:
    """Normalize voucher to NCBI 'institution:specimen-id' format."""
    v = str(voucher).strip()
    if not v or v == "nan" or v.lower() == "no voucher":
        return ""

    # Remove type designations
    v = re.sub(r"\s+(Holotype|Paratype)$", "", v, flags=re.IGNORECASE)

    # Already has colon separator
    if ":" in v:
        return v

    # "SMNH 211581" → "SMNH:211581"
    parts = v.split(None, 1)
    if len(parts) == 2 and parts[0].isalpha():
        return f"{parts[0]}:{parts[1]}"

    # "RBINS IG 33953/09.132.01" → "RBINS IG:33953/09.132.01"
    m = re.match(r"^([A-Z]+)\s+([A-Z]+)\s+(.+)$", v)
    if m:
        return f"{m.group(1)} {m.group(2)}:{m.group(3)}"

    # "NSMT An-1915-1917" → "NSMT:An-1915-1917"
    if parts[0].isalpha():
        return f"{parts[0]}:{parts[1]}"

    return v


def load_literature() -> dict:
    """Load literature Excel into accession → metadata lookup."""
    df = pd.read_excel(EXCEL, sheet_name="Blad1", header=2)
    df.columns = [
        "Species", "Specimen_no", "Voucher", "Country", "State",
        "Municipality", "Locality", "Habitat", "N", "E",
        "Coll_Date", "Leg", "12S", "16S", "18S", "28S", "COI",
    ]
    df = df.iloc[1:].reset_index(drop=True)

    # Build specimen-level metadata
    specimens = {}
    for _, row in df.iterrows():
        meta = {
            "country": str(row["Country"]).strip(),
            "state": str(row["State"]).strip(),
            "municipality": str(row["Municipality"]).strip(),
            "locality": str(row["Locality"]).strip(),
            "habitat": str(row["Habitat"]).strip(),
            "n": row["N"],
            "e": row["E"],
            "coll_date": str(row["Coll_Date"]).strip(),
            "leg": str(row["Leg"]).strip(),
            "voucher": str(row["Voucher"]).strip(),
            "species": str(row["Species"]).strip(),
        }
        # Each specimen maps to multiple accessions via gene columns
        for gene in ["12S", "16S", "18S", "28S", "COI"]:
            acc = str(row[gene]).strip()
            if acc and acc != "nan" and acc != "-":
                specimens[acc] = meta

    return specimens


def load_gb_state() -> dict:
    """Load current GB metadata state from g2t final.csv."""
    final = pd.read_csv(FINAL_CSV)
    lookup = {}
    for _, row in final.iterrows():
        acc = str(row["ACCESSION"]).strip()
        lookup[acc] = {
            "lat_lon": row.get("lat_lon"),
            "collection_date": row.get("collection_date"),
            "collected_by": row.get("collected_by"),
            "geo_loc_name": row.get("geo_loc_name"),
            "specimen_voucher": row.get("specimen_voucher"),
            "isolate": row.get("isolate"),
            "isolation_source": row.get("isolation_source", ""),
        }
    return lookup


def main():
    lit = load_literature()
    gb = load_gb_state()

    print(f"文献 accessions: {len(lit)}")
    print(f"GB final.csv: {len(gb)} accessions")

    # Build update rows
    ncbi_columns = [
        "acc. num.",
        "geo_loc_name",
        "Lat_Lon",
        "Collection_date",
        "Collected_by",
        "Specimen_voucher",
        "Isolation_source",
    ]

    rows = []
    stats = {"matched": 0, "updated": 0, "fields_filled": {c: 0 for c in ncbi_columns[1:]}}

    for acc, lit_meta in sorted(lit.items()):
        if acc not in gb:
            continue
        stats["matched"] += 1
        gb_meta = gb[acc]

        row = {"acc. num.": acc}
        any_update = False

        # Country (geo_loc_name) — only if GB lacks it
        if is_empty(gb_meta["geo_loc_name"]):
            val = format_geo_loc_name(
                lit_meta["country"], lit_meta["state"],
                lit_meta["municipality"], lit_meta["locality"],
            )
            if val:
                row["Country (geo_loc_name)"] = val
                any_update = True
                stats["fields_filled"]["geo_loc_name"] += 1

        # Lat_Lon
        if is_empty(gb_meta["lat_lon"]):
            val = format_lat_lon(lit_meta["n"], lit_meta["e"])
            if val:
                row["Lat_Lon"] = val
                any_update = True
                stats["fields_filled"]["Lat_Lon"] += 1

        # Collection_date
        if is_empty(gb_meta["collection_date"]):
            val = normalize_date(lit_meta["coll_date"])
            if val:
                row["Collection_date"] = val
                any_update = True
                stats["fields_filled"]["Collection_date"] += 1

        # Collected_by
        if is_empty(gb_meta["collected_by"]):
            val = lit_meta["leg"]
            if val and val != "nan":
                row["Collected_by"] = val
                any_update = True
                stats["fields_filled"]["Collected_by"] += 1

        # Specimen_voucher
        if is_empty(gb_meta["specimen_voucher"]) and is_empty(gb_meta["isolate"]):
            val = format_voucher(lit_meta["voucher"])
            if val:
                row["Specimen_voucher"] = val
                any_update = True
                stats["fields_filled"]["Specimen_voucher"] += 1

        # Isolation_source (habitat)
        if is_empty(gb_meta.get("isolation_source")):
            val = lit_meta["habitat"]
            if val and val != "nan":
                row["Isolation_source"] = val
                any_update = True
                stats["fields_filled"]["Isolation_source"] += 1

        if any_update:
            rows.append(row)
            stats["updated"] += 1

    # Write TSV
    out_df = pd.DataFrame(rows, columns=ncbi_columns)
    out_df = out_df.fillna("")
    out_df.to_csv(OUTPUT, sep="\t", index=False)

    # Report
    print(f"\n{'=' * 60}")
    print(f"NCBI Source Modifiers 更新表")
    print(f"{'=' * 60}")
    print(f"文献 → GB 匹配: {stats['matched']} accessions")
    print(f"需要更新: {stats['updated']} accessions")
    print()
    print(f"各字段补充数量:")
    for col in ncbi_columns[1:]:
        print(f"  {col:30s} {stats['fields_filled'][col]:5d} 条")
    print()
    print(f"输出: {OUTPUT}")
    print(f"格式: Tab-delimited, NCBI Source Modifiers Table")

    # Preview first 3 rows
    if len(out_df) > 0:
        print(f"\n前 3 行预览:")
        for _, r in out_df.head(3).iterrows():
            print(f"  {r['acc. num.']}")
            for col in ncbi_columns[1:]:
                val = r[col]
                if val:
                    print(f"    {col}: {val}")


if __name__ == "__main__":
    main()
