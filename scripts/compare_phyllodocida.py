#!/usr/bin/env python3
"""
Compare g2t Phyllodocida output with reference data from literature.

Usage:
    python scripts/compare_phyllodocida.py
    python scripts/compare_phyllodocida.py --no-mtgenome-fill
    python scripts/compare_phyllodocida.py --ncbi-lookup

Compares:
    - g2t organized_species_voucher.csv (with mtgenome propagation)
    - Reference: ceshi2/TableS1_Oct2024.xlsx (88 species x 7 gene types)
"""

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd

BASE = Path("/media/deyuan/217ce44c-b5de-45ed-8720-deebdff85ece/CLAUDE/ceshi")
G2T_OUTPUT = Path("Phyllodocida_output")

GENE_MAP = {
    "COX1": "coi",
    "16S": "16s",
    "18S": "18s",
    "28S": "28s",
    "H3": "h3",
    "CYTB": "cob",
    "Mitogenome": "mtgenome",
}

# Genus remapping: literature name -> {species_epithet: new_genus}
# Based on Zhang et al. 2024 and Hatch & Rouse 2020 reclassifications
GENUS_REMAP = {
    "branchinotogluma": {
        "sandersi": "Cladopolynoe",
        "tunnicliffeae": "Cladopolynoe",
        "jiaolongae": "Cladopolynoe",
        "trifurcus": "Cladopolynoe",
        "marianus": "Branchipolynoe",
        "segonzaci": "Branchipolynoe",
        "kaireiensis": "Photinopolynoe",
        "robusta": "Photinopolynoe",
        "nanhaiensis": "Photinopolynoe",
        "bipapillata": "Photinopolynoe",
        "elytropapillata": "Photinopolynoe",
        "hessleri": "Photinopolynoe",
        "japonicus": "Photinopolynoe",
        "nikkoensis": "Photinopolynoe",
        "ovata": "Photinopolynoe",
        "pettiboneae": "Photinopolynoe",
        "sagamiensis": "Photinopolynoe",
    },
    "lepidonotopodium": {
        "piscesae": "Mamiwata",
        "williamsae": "Mamiwata",
        "riftense": "Mamiwata",
    },
    "levensteuniella": {
        "intermedia": "Themis",
        "longqiensis": "Themis",
        "manusensis": "Themis",
        "pettiboneae": "Themis",
        "plicata": "Themis",
        "undomarginata": "Themis",
    },
}


def normalize_species(name: str) -> str:
    name = name.strip()
    # Remove cf., sp. nov., sensu, gen. nov.
    name = re.sub(r"\b(cf\.|sp\.?\s*nov\.?\s*\d*|gen\.?\s*nov\.?|sensu\s+\S+)\b", " ", name, flags=re.I)
    # Remove (Author, Year) parentheses
    name = re.sub(r"\([^)]*\d{4}[^)]*\)", "", name)
    # Remove trailing author citations: "Han, 2023" or "Han, Zhou, Chen & Wang, 2023"
    name = re.sub(r",?\s*[A-Z][a-z]+.*&.*\d{4}.*$", "", name)
    name = re.sub(r",?\s*[A-Z][a-z]+,?\s*\d{4}.*$", "", name)
    name = re.sub(r"\s+", " ", name).strip()
    return name.lower()


def extract_genus_species(name: str) -> tuple[str, str]:
    norm = normalize_species(name)
    parts = norm.split()
    return (parts[0] if parts else ""), (parts[1] if len(parts) > 1 else "")


def build_remap_candidates(lit_name: str) -> list[str]:
    """Generate candidate names: original + remapped genus versions."""
    genus, species = extract_genus_species(lit_name)
    candidates = [lit_name]

    genus_lower = genus.lower()
    if genus_lower not in GENUS_REMAP:
        return candidates

    species_map = GENUS_REMAP[genus_lower]
    new_genera = set()

    if species and species in species_map:
        new_genera.add(species_map[species])
    else:
        # For sp. nov. / cf. names, try all possible new genera
        new_genera = set(species_map.values())

    for new_genus in new_genera:
        if species:
            candidates.append(f"{new_genus} {species}")

    return candidates


def match_species(lit_species: str, g2t_organisms: dict[str, str]) -> tuple[str | None, str]:
    """Match literature species to g2t organism via normalized name comparison.

    g2t_organisms maps normalized_name -> original_name.
    Returns (original_organism, method).
    """
    candidates = build_remap_candidates(lit_species)
    is_remap = False

    for candidate in candidates:
        cand_norm = normalize_species(candidate)
        cand_parts = cand_norm.split()

        # Exact match
        if cand_norm in g2t_organisms:
            method = "remap" if is_remap else "exact"
            return g2t_organisms[cand_norm], method

        # Genus + species prefix match (first 4 chars of epithet)
        # Skip if candidate is "cf." or "sp." — too ambiguous for prefix match
        if len(cand_parts) >= 2 and cand_parts[1] not in ("cf", "sp"):
            prefix = f"{cand_parts[0]} {cand_parts[1][:4]}"
            for norm_org in g2t_organisms:
                org_parts = norm_org.split()
                if len(org_parts) < 2 or org_parts[1] in ("cf", "sp"):
                    continue
                if f"{org_parts[0]} {org_parts[1][:4]}" == prefix:
                    method = "remap" if is_remap else "prefix"
                    return g2t_organisms[norm_org], method

        is_remap = True

    return None, "none"


def accession_match(lit_acc: str, g2t_accs: str) -> bool:
    if not lit_acc or not g2t_accs:
        return False
    lit_acc = lit_acc.strip()
    for g2t_acc in g2t_accs.split("; "):
        g2t_acc = g2t_acc.strip()
        if g2t_acc == lit_acc or g2t_acc.startswith(lit_acc + "."):
            return True
    return False


def build_accession_organism_map(assigned_path: Path) -> dict[str, str]:
    """Build accession -> organism lookup from assigned_genes_types_all.csv."""
    if not assigned_path.exists():
        return {}
    df = pd.read_csv(assigned_path, low_memory=False, usecols=["ACCESSION", "organism"])
    acc_map = {}
    for _, r in df.dropna(subset=["ACCESSION", "organism"]).iterrows():
        acc = str(r["ACCESSION"]).strip().upper()
        acc_map[acc] = str(r["organism"])
    return acc_map


def ncbi_lookup_accession(accession: str, email: str, api_key: str = None) -> dict:
    try:
        from Bio import Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        text = handle.read()
        handle.close()

        organism = ""
        for line in text.split("\n"):
            if line.startswith("  ORGANISM"):
                organism = line.replace("  ORGANISM", "").strip()
                break
        if not organism:
            for line in text.split("\n"):
                if line.startswith("SOURCE"):
                    organism = line.replace("SOURCE", "").strip()
                    break
        return {"accession": accession, "organism": organism}
    except Exception as e:
        return {"accession": accession, "organism": "", "error": str(e)}


def parse_args():
    p = argparse.ArgumentParser(description="Compare g2t output with literature reference")
    p.add_argument("--no-mtgenome-fill", action="store_true",
                   help="Disable mtgenome fill for missing sub-genes")
    p.add_argument("--ncbi-lookup", action="store_true",
                   help="Do NCBI reverse lookup for still-unmatched species")
    p.add_argument("--email", default="deyuanyang@example.com")
    p.add_argument("--api-key", default="62733a242c3b7a3b0b93e6d65e108798b708")
    return p.parse_args()


def main():
    args = parse_args()
    mtgenome_fill = not args.no_mtgenome_fill

    ref_path = BASE / "ceshi2" / "TableS1_Oct2024.xlsx"
    organized_path = G2T_OUTPUT / "organized_genes" / "organized_species_voucher.csv"
    assigned_path = G2T_OUTPUT / "labeled_genes" / "assigned_genes_types_all.csv"

    if not organized_path.exists():
        print("Error: g2t organized output not found. Run g2t first.")
        sys.exit(1)

    ref_df = pd.read_excel(ref_path)
    print(f"Reference: {len(ref_df)} species from TableS1_Oct2024.xlsx")

    org_df = pd.read_csv(organized_path, low_memory=False)
    print(f"g2t organized: {len(org_df)} species vouchers")

    gene_cols = ["mtgenome", "coi", "16s", "12s", "cob", "cox2", "cox3", "18s", "28s", "h3"]
    for col in gene_cols:
        if col in org_df.columns:
            org_df[col] = org_df[col].fillna("").astype(str).replace({"nan": ""})

    # Aggregate gene data per organism (combine all vouchers)
    org_gene_data = {}
    for _, r in org_df.iterrows():
        org = r["organism"]
        if org not in org_gene_data:
            org_gene_data[org] = {col: set() for col in gene_cols}
        for col in gene_cols:
            val = str(r.get(col, ""))
            if val and val != "nan":
                for acc in val.split("; "):
                    acc = acc.strip()
                    if acc:
                        org_gene_data[org][col].add(acc)

    # Convert sets to sorted strings
    org_lookup = {}
    for org, genes in org_gene_data.items():
        row_data = {}
        for col in gene_cols:
            row_data[col] = "; ".join(sorted(genes[col]))
        org_lookup[org] = row_data

    # Normalized organism index for fast matching
    g2t_norm_map = {}
    for org in org_lookup:
        norm = normalize_species(org)
        if norm not in g2t_norm_map:
            g2t_norm_map[norm] = org

    print(f"g2t unique organisms: {len(g2t_norm_map)}")

    # Build accession -> organism map for local reverse lookup
    print("Building accession index...")
    acc_org_map = build_accession_organism_map(assigned_path)
    print(f"Accession index: {len(acc_org_map)} entries")

    gene_labels = [g for g in GENE_MAP if g in ref_df.columns]

    # ---------- Pass 1: Name-based matching + gene comparison ----------
    results = []
    total_lit = 0
    total_found = 0
    total_acc_match = 0
    total_propagated = 0

    for _, row in ref_df.iterrows():
        sp = row["Species"]
        matched_org, match_method = match_species(sp, g2t_norm_map)
        org_data = org_lookup.get(matched_org) if matched_org else None

        record = {
            "species": sp,
            "g2t_organism": matched_org or "NOT MATCHED",
            "match_method": match_method,
        }

        has_any_lit = False
        all_match = True
        lit_count = 0
        found_count = 0
        acc_match_count = 0
        propagated_count = 0
        gene_covered = 0  # g2t has any data for this gene (even if different accession)

        for glit in gene_labels:
            g2t_key = GENE_MAP[glit]
            lit_acc = str(row.get(glit, "")).strip().strip("﻿﻿")
            if lit_acc in ("-", "nan", ""):
                lit_acc = ""

            g2t_accs = ""
            if org_data is not None:
                g2t_accs = org_data.get(g2t_key, "")

            was_propagated = False
            if mtgenome_fill and g2t_key != "mtgenome" and not g2t_accs and org_data is not None:
                mtg_accs = org_data.get("mtgenome", "")
                if mtg_accs:
                    g2t_accs = mtg_accs
                    was_propagated = True

            record[f"lit_{glit}"] = lit_acc
            record[f"g2t_{g2t_key}"] = g2t_accs[:200] if g2t_accs else ""

            if lit_acc:
                has_any_lit = True
                lit_count += 1
                if g2t_accs:
                    found_count += 1
                    gene_covered += 1
                if accession_match(lit_acc, g2t_accs):
                    acc_match_count += 1
                    total_acc_match += 1
                else:
                    all_match = False
                    if was_propagated:
                        propagated_count += 1
            elif g2t_accs:
                found_count += 1

        total_lit += lit_count
        total_found += found_count
        total_propagated += propagated_count

        record["lit_genes"] = lit_count
        record["g2t_found"] = found_count
        record["acc_matched"] = acc_match_count
        record["gene_covered"] = gene_covered
        record["propagated"] = propagated_count

        if not has_any_lit:
            record["match_type"] = "NO_LIT_DATA"
        elif matched_org is None:
            record["match_type"] = "NOT_FOUND"
        elif all_match and acc_match_count == lit_count:
            record["match_type"] = "MATCH"
        else:
            record["match_type"] = "PARTIAL"

        if match_method == "remap":
            record["detail"] = f"Remapped: {sp.split()[0]} -> {matched_org.split()[0]}"
        elif matched_org is None:
            record["detail"] = "Species not found in g2t data"
        elif record["match_type"] == "MATCH":
            record["detail"] = "All accessions matched"
        else:
            record["detail"] = f"{acc_match_count}/{lit_count} accessions matched"

        results.append(record)

    result_df = pd.DataFrame(results)

    # ---------- Pass 2: Local accession-based lookup for NOT_FOUND ----------
    still_not_found = result_df[result_df["match_type"] == "NOT_FOUND"]
    if len(still_not_found) > 0:
        print(f"\n--- Local Accession Lookup for {len(still_not_found)} NOT_FOUND species ---")
        for idx in still_not_found.index:
            r = result_df.loc[idx]
            sp = r["species"]
            found_org = None
            for glit in gene_labels:
                lit_acc = str(r.get(f"lit_{glit}", "")).strip()
                if not lit_acc or lit_acc in ("", "nan", "-"):
                    continue
                # Strip asterisk (pre-publication marker)
                clean_acc = lit_acc.rstrip("*").strip()
                # Try with and without version suffix
                for lookup_acc in [clean_acc, clean_acc.split(".")[0]]:
                    if lookup_acc.upper() in acc_org_map:
                        found_org = acc_org_map[lookup_acc.upper()]
                        break
                    # Also try with .1 suffix
                    if f"{lookup_acc.upper()}.1" in acc_org_map:
                        found_org = acc_org_map[f"{lookup_acc.upper()}.1"]
                        break
                if found_org:
                    break

            if found_org:
                print(f"  {sp}")
                print(f"    -> Found via accession as: {found_org}")
                result_df.at[idx, "g2t_organism"] = found_org
                result_df.at[idx, "match_method"] = "accession_lookup"
                result_df.at[idx, "match_type"] = "ACCESSION_RESOLVED"
                result_df.at[idx, "detail"] = f"Found via accession lookup: {found_org}"

                # Re-do gene comparison with resolved organism
                org_data = org_lookup.get(found_org)
                if org_data is not None:
                    acc_match_count = 0
                    lit_count = 0
                    for glit in gene_labels:
                        g2t_key = GENE_MAP[glit]
                        lit_acc = str(ref_df.loc[ref_df["Species"] == sp].iloc[0].get(glit, "")).strip().strip("﻿﻿")
                        if lit_acc in ("-", "nan", ""):
                            lit_acc = ""
                        lit_acc = lit_acc.rstrip("*").strip()
                        g2t_accs = org_data.get(g2t_key, "")

                        if mtgenome_fill and g2t_key != "mtgenome" and not g2t_accs:
                            mtg_accs = org_data.get("mtgenome", "")
                            if mtg_accs:
                                g2t_accs = mtg_accs

                        if lit_acc:
                            lit_count += 1
                            if accession_match(lit_acc, g2t_accs):
                                acc_match_count += 1
                        result_df.at[idx, f"g2t_{g2t_key}"] = g2t_accs[:200]

                    result_df.at[idx, "acc_matched"] = acc_match_count
                    result_df.at[idx, "lit_genes"] = lit_count
                    if acc_match_count == lit_count and lit_count > 0:
                        result_df.at[idx, "match_type"] = "MATCH"
                    elif acc_match_count > 0:
                        result_df.at[idx, "match_type"] = "PARTIAL"
            else:
                # Check if accessions have asterisks (pre-publication)
                has_asterisk = False
                for glit in gene_labels:
                    v = str(r.get(f"lit_{glit}", ""))
                    if "*" in v:
                        has_asterisk = True
                        break
                if has_asterisk:
                    result_df.at[idx, "detail"] = "Pre-publication accessions (not yet in GenBank)"
                else:
                    result_df.at[idx, "detail"] = "Not found locally or by accession"

    # ---------- Pass 3: NCBI reverse lookup (optional) ----------
    if args.ncbi_lookup:
        remaining = result_df[result_df["match_type"] == "NOT_FOUND"]
        if len(remaining) > 0:
            print(f"\n--- NCBI Reverse Lookup for {len(remaining)} remaining species ---")
            for idx in remaining.index:
                r = result_df.loc[idx]
                sp = r["species"]
                for glit in gene_labels:
                    lit_acc = str(r.get(f"lit_{glit}", "")).strip().rstrip("*")
                    if not lit_acc or lit_acc in ("", "nan", "-"):
                        continue
                    print(f"  Looking up {lit_acc} for {sp}...")
                    info = ncbi_lookup_accession(lit_acc, args.email, args.api_key)
                    if info and info.get("organism"):
                        ncbi_org = info["organism"]
                        print(f"    NCBI: {lit_acc} -> {ncbi_org}")
                        ncbi_matched, _ = match_species(ncbi_org, g2t_norm_map)
                        if ncbi_matched:
                            print(f"    -> Found in g2t as: {ncbi_matched}")
                            result_df.at[idx, "g2t_organism"] = ncbi_matched
                            result_df.at[idx, "match_type"] = "NCBI_RESOLVED"
                            result_df.at[idx, "detail"] = f"NCBI: {ncbi_org} -> {ncbi_matched}"
                        else:
                            print(f"    -> Not in g2t: {ncbi_org}")
                    time.sleep(0.15)
                    break

    # ---------- Output ----------
    out_path = G2T_OUTPUT / "comparison_results.csv"
    result_df.to_csv(out_path, index=False)

    # ---------- Build detailed line-by-line report ----------
    report_lines = []
    total_ok = 0
    total_miss = 0
    total_extra = 0  # g2t has data but lit doesn't

    for _, r in result_df.iterrows():
        sp = r["species"]
        org = r["g2t_organism"]
        mm = r["match_method"]

        # Species header
        name_note = ""
        if org != "NOT MATCHED" and normalize_species(sp) != normalize_species(org):
            name_note = f"  (g2t name: {org})"
        report_lines.append(f"")
        report_lines.append(f"{'='*70}")
        report_lines.append(f"{sp}{name_note}")
        report_lines.append(f"{'='*70}")

        for glit in gene_labels:
            g2t_key = GENE_MAP[glit]
            lit_acc = str(r.get(f"lit_{glit}", "")).strip().strip("﻿﻿")
            if lit_acc in ("-", "nan", ""):
                lit_acc = ""
            g2t_accs = str(r.get(f"g2t_{g2t_key}", "")).strip()
            if g2t_accs in ("-", "nan", ""):
                g2t_accs = ""

            if not lit_acc and not g2t_accs:
                continue  # neither side has data

            if lit_acc and accession_match(lit_acc, g2t_accs):
                total_ok += 1
                report_lines.append(f"  {glit:12s}  OK   lit={lit_acc}")
            elif lit_acc:
                total_miss += 1
                # Determine reason
                if "*" in lit_acc:
                    reason = "pre-publication (not yet in GenBank)"
                elif not g2t_accs:
                    reason = "g2t has no data for this gene"
                else:
                    reason = "g2t has different accession (same gene, different sequence)"
                lit_display = lit_acc
                g2t_display = g2t_accs[:60] if g2t_accs else "(none)"
                report_lines.append(f"  {glit:12s}  MISS lit={lit_display}")
                report_lines.append(f"  {'':12s}       g2t={g2t_display}")
                report_lines.append(f"  {'':12s}       reason: {reason}")
            else:
                # lit has no data but g2t does
                total_extra += 1
                report_lines.append(f"  {glit:12s}  EXTRA (g2t only)  g2t={g2t_accs[:50]}")

    # Summary header
    summary = []
    summary.append(f"")
    summary.append(f"{'#'*70}")
    summary.append(f"# Comparison Report: 88 reference species vs g2t")
    summary.append(f"{'#'*70}")
    summary.append(f"")
    summary.append(f"Legend:")
    summary.append(f"  OK   = literature accession found in g2t")
    summary.append(f"  MISS = literature accession NOT found in g2t")
    summary.append(f"  EXTRA = g2t has data but literature does not")
    summary.append(f"")
    summary.append(f"Totals:  {total_ok} OK,  {total_miss} MISS,  {total_extra} EXTRA")
    summary.append(f"Match rate: {total_ok}/{total_ok+total_miss} ({total_ok/(total_ok+total_miss)*100:.1f}%)")
    summary.append(f"")

    # Print to terminal
    for line in summary:
        print(line)
    for line in report_lines:
        print(line)

    # Save report to file
    report_path = G2T_OUTPUT / "comparison_report.txt"
    with open(report_path, "w") as f:
        f.write("\n".join(summary + report_lines))
    print(f"\nReport saved: {report_path}")
    print(f"CSV saved:    {out_path}")


if __name__ == "__main__":
    main()
