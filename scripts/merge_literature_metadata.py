#!/usr/bin/env python3
"""
Merge literature metadata into g2t final.csv output.

Supplements GB-extracted metadata with richer literature data (GPS, dates, collectors, etc.)
GB values take priority — literature only fills empty fields.

Usage:
    python scripts/merge_literature_metadata.py
"""

import pandas as pd
from pathlib import Path

BASE = Path("/media/deyuan/217ce44c-b5de-45ed-8720-deebdff85ece/CLAUDE/ceshi/ceshi5")
EXCEL = BASE / "zlae141_suppl_supplementary_tables_s1.xlsx"
FINAL_CSV = BASE / "haplotaxis_output3/gb_metadata/final.csv"
OUTPUT = BASE / "haplotaxis_output3/gb_metadata/final_enriched.csv"


def load_literature() -> pd.DataFrame:
    df = pd.read_excel(EXCEL, sheet_name="Blad1", header=2)
    df.columns = [
        "Species", "Specimen_no", "Voucher", "Country", "State",
        "Municipality", "Locality", "Habitat", "N", "E",
        "Coll_Date", "Leg", "12S", "16S", "18S", "28S", "COI",
    ]
    df = df.iloc[1:].reset_index(drop=True)

    # Melt: one row per accession-gene pair
    records = []
    for _, row in df.iterrows():
        base = {
            "lit_species": str(row["Species"]).strip(),
            "lit_specimen_no": str(row["Specimen_no"]).strip(),
            "lit_voucher": str(row["Voucher"]).strip(),
            "lit_country": str(row["Country"]).strip(),
            "lit_state": str(row["State"]).strip(),
            "lit_municipality": str(row["Municipality"]).strip(),
            "lit_locality": str(row["Locality"]).strip(),
            "lit_habitat": str(row["Habitat"]).strip(),
            "lit_lat": row["N"],
            "lit_lon": row["E"],
            "lit_coll_date": str(row["Coll_Date"]).strip(),
            "lit_leg": str(row["Leg"]).strip(),
        }
        for gene in ["12S", "16S", "18S", "28S", "COI"]:
            acc = str(row[gene]).strip()
            if acc and acc != "nan":
                rec = dict(base, lit_gene=gene)
                records.append((acc, rec))

    # Build lookup: accession -> literature data
    lookup = {}
    for acc, rec in records:
        if acc not in lookup:
            lookup[acc] = rec
    return lookup


def main():
    final = pd.read_csv(FINAL_CSV)
    lit_lookup = load_literature()

    print(f"g2t final.csv: {len(final)} records")
    print(f"Literature lookup: {len(lit_lookup)} accessions")

    # New columns to add
    lit_cols = [
        "lit_country", "lit_state", "lit_municipality", "lit_locality",
        "lit_habitat", "lit_lat", "lit_lon", "lit_coll_date", "lit_leg",
        "lit_voucher", "lit_species",
    ]
    for col in lit_cols:
        final[col] = ""

    matched = 0
    for idx, row in final.iterrows():
        acc = str(row["ACCESSION"]).strip()
        if acc in lit_lookup:
            matched += 1
            lit = lit_lookup[acc]
            for col in lit_cols:
                final.at[idx, col] = lit.get(col, "")

    # Fill GPS from literature where lat_lon is empty
    filled_gps = 0
    filled_date = 0
    filled_collector = 0
    filled_habitat = 0
    for idx, row in final.iterrows():
        # GPS
        if (not row.get("lat_lon") or str(row["lat_lon"]) in ("", "nan")):
            if row.get("lit_lat") and str(row["lit_lat"]) not in ("", "nan"):
                lat = row["lit_lat"]
                lon = row["lit_lon"]
                final.at[idx, "lat_lon"] = f"{lat} {lon}"
                filled_gps += 1
        # Collection date
        if (not row.get("collection_date") or str(row["collection_date"]) in ("", "nan")):
            if row.get("lit_coll_date") and str(row["lit_coll_date"]) not in ("", "nan"):
                final.at[idx, "collection_date"] = row["lit_coll_date"]
                filled_date += 1
        # Collector
        if (not row.get("collected_by") or str(row["collected_by"]) in ("", "nan")):
            if row.get("lit_leg") and str(row["lit_leg"]) not in ("", "nan"):
                final.at[idx, "collected_by"] = row["lit_leg"]
                filled_collector += 1

    final.to_csv(OUTPUT, index=False)

    # Report
    print(f"\n{'=' * 60}")
    print(f"合并结果")
    print(f"{'=' * 60}")
    print(f"匹配: {matched}/{len(final)} accessions 在文献中找到")
    print()
    print(f"字段补充:")
    print(f"  GPS 坐标 (lat_lon):       补充了 {filled_gps} 条")
    print(f"  采集日期 (collection_date): 补充了 {filled_date} 条")
    print(f"  采集人 (collected_by):     补充了 {filled_collector} 条")
    print()

    # Before/after comparison
    print(f"{'字段':25s} {'合并前':>10s} {'合并后':>10s} {'提升':>10s}")
    print("-" * 60)

    def count_nonempty(df, col):
        return (df[col].notna() & (df[col].astype(str) != "") & (df[col].astype(str) != "nan")).sum()

    # Recount after merge
    final2 = pd.read_csv(OUTPUT)
    comparisons = [
        ("lat_lon", "lat_lon"),
        ("collection_date", "collection_date"),
        ("collected_by", "collected_by"),
        ("lit_country", "lit_country"),
        ("lit_state", "lit_state"),
        ("lit_municipality", "lit_municipality"),
        ("lit_locality", "lit_locality"),
        ("lit_habitat", "lit_habitat"),
    ]
    # GB original counts (re-read original)
    orig = pd.read_csv(FINAL_CSV)
    for label, col in comparisons:
        before = count_nonempty(orig, col) if col in orig.columns else 0
        after = count_nonempty(final2, col)
        delta = after - before
        print(f"  {label:23s} {before:8d}   {after:8d}   {delta:+8d}")

    print(f"\n输出: {OUTPUT}")
    print(f"列数: {len(final2.columns)} (原 {len(orig.columns)})")


if __name__ == "__main__":
    main()
