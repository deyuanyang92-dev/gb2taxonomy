#!/usr/bin/env python3
"""
Build g2t v2.0 gene dictionary by merging SynGenes TSV data with existing g2t synonyms.

Usage:
    python scripts/build_gene_dict_v2.py [--output src/g2t/data/gene_dict_v2.yaml]

Reads:
    - src/g2t/data/gene_dict.yaml (existing g2t v1 dictionary)
    - ../SynGenes/data Names/mt/*.tsv (18 mitochondrial gene files)
    - ../SynGenes/data Names/cp/*.tsv (101 chloroplast gene files)

Outputs:
    - gene_dict_v2.yaml (hierarchical YAML with category/scope fields)
"""

import os
import sys
import argparse
import csv
from pathlib import Path
from collections import defaultdict

# Try to import yaml, fall back to manual output
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


# ============================================================
# Gene name mappings between SynGenes and g2t naming conventions
# ============================================================

# SynGenes Short Name -> g2t gene key (for overlapping genes)
SYNGENES_TO_G2T_MAP = {
    "COI": "coi",
    "COII": "cox2",
    "COIII": "cox3",
    "CYTB": "cob",
    "12S": "12s",
    "16S": "16s",
}

# SynGenes Short Name -> g2t gene key (new mitochondrial genes)
NEW_MT_GENES = {
    "ND1": "nd1",
    "ND2": "nd2",
    "ND3": "nd3",
    "ND4": "nd4",
    "ND4L": "nd4l",
    "ND5": "nd5",
    "ND6": "nd6",
    "ATP6": "atp6",
    "ATP8": "atp8",
    "Control Region": "control_region",
    "OH": "oh",
    "OL": "ol",
}

# All mitochondrial mapping combined
MT_MAP = {**SYNGENES_TO_G2T_MAP, **NEW_MT_GENES}

# Synonyms to skip because they're too short/ambiguous and cause false positives
SKIP_SYNONYMS = {
    "I", "COX", "CO", "ND", "ATP", "rbcL",  # too ambiguous alone
}

# Minimum synonym length to include (shorter ones cause false positives)
MIN_SYNONYM_LEN = 3


def find_g2t_root():
    """Find g2t project root."""
    script_dir = Path(__file__).resolve().parent
    return script_dir.parent


def find_syngenes_root():
    """Find SynGenes data root."""
    g2t_root = find_g2t_root()
    return g2t_root.parent / "SynGenes"


def read_g2t_gene_dict(yaml_path):
    """Read existing g2t v1 gene dictionary."""
    if not os.path.exists(yaml_path):
        print(f"Warning: {yaml_path} not found, starting fresh")
        return {}

    if HAS_YAML:
        with open(yaml_path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
        result = {}
        for k, v in data.get('gene_types', {}).items():
            synonyms = v.get('synonyms', []) if isinstance(v, dict) else (v if isinstance(v, list) else [])
            result[k.lower()] = [str(s) for s in synonyms if isinstance(s, str)]
        return result
    else:
        print("Error: PyYAML required. Install with: pip install pyyaml")
        sys.exit(1)


def read_syngenes_tsv(tsv_path):
    """Read a SynGenes TSV file and return list of Full Name synonyms."""
    synonyms = []
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)  # skip header
        for row in reader:
            if len(row) >= 1 and row[0].strip():
                name = row[0].strip()
                if name and name not in SKIP_SYNONYMS and len(name) >= MIN_SYNONYM_LEN:
                    synonyms.append(name)
    return synonyms


def classify_existing_genes(gene_key):
    """Determine category and scope for existing g2t v1 genes."""
    mito_genes = {"coi", "16s", "12s", "cob", "mtgenome", "cox2", "cox3"}
    nuclear_genes = {"18-28s", "its1-its2", "ef-1", "28s", "18s", "h3"}

    if gene_key == "mtgenome":
        return "mitochondrial", ["metazoa", "fungi"]
    elif gene_key in mito_genes:
        return "mitochondrial", ["metazoa", "fungi"]
    elif gene_key in nuclear_genes:
        return "nuclear", ["metazoa", "viridiplantae", "fungi"]
    else:
        return "unknown", ["all"]


def build_v2_dict(g2t_dict, mt_dir, cp_dir):
    """Build the v2 gene dictionary with category/scope metadata."""
    v2 = {}

    # === Existing g2t genes (enriched with SynGenes synonyms) ===
    for gene_key, existing_synonyms in g2t_dict.items():
        category, scope = classify_existing_genes(gene_key)

        # Check if SynGenes has additional synonyms for this gene
        syngenes_synonyms = []
        for syn_short, g2t_key in MT_MAP.items():
            if g2t_key == gene_key:
                tsv_file = mt_dir / f"{syn_short}.tsv"
                if tsv_file.exists():
                    syngenes_synonyms = read_syngenes_tsv(tsv_file)
                break

        # Merge synonyms (deduplicate case-insensitively, preserve original case)
        seen_lower = set()
        merged = []
        for s in existing_synonyms + syngenes_synonyms:
            s_lower = s.lower()
            if s_lower not in seen_lower and s not in SKIP_SYNONYMS and len(s) >= MIN_SYNONYM_LEN:
                seen_lower.add(s_lower)
                merged.append(s)

        v2[gene_key] = {
            "category": category,
            "scope": scope,
            "synonyms": merged,
        }

    # === New mitochondrial genes from SynGenes ===
    for syn_short, g2t_key in NEW_MT_GENES.items():
        if g2t_key in v2:
            continue  # already processed above
        tsv_file = mt_dir / f"{syn_short}.tsv"
        if tsv_file.exists():
            synonyms = read_syngenes_tsv(tsv_file)
            v2[g2t_key] = {
                "category": "mitochondrial",
                "scope": ["metazoa", "fungi"],
                "synonyms": synonyms,
            }

    # === Chloroplast genes from SynGenes ===
    if cp_dir.exists():
        for tsv_file in sorted(cp_dir.glob("*.tsv")):
            gene_name = tsv_file.stem  # e.g., "rbcL", "matK", "psbA"
            g2t_key = gene_name.lower()

            # Avoid collision with existing genes
            if g2t_key in v2:
                g2t_key = f"cp_{g2t_key}"

            synonyms = read_syngenes_tsv(tsv_file)
            if synonyms:
                v2[g2t_key] = {
                    "category": "chloroplast",
                    "scope": ["viridiplantae"],
                    "synonyms": synonyms,
                }

    return v2


def write_v2_yaml(v2_dict, output_path):
    """Write the v2 dictionary to YAML."""
    # Build output structure
    output = {
        "version": "2.0",
        "metadata": {
            "source": "g2t v2.0 + SynGenes integration",
            "description": "Unified gene dictionary with organism-aware scope",
        },
        "gene_types": {},
    }

    # Sort: mitochondrial first, then chloroplast, then nuclear, then unknown
    category_order = {"mitochondrial": 0, "chloroplast": 1, "nuclear": 2, "unknown": 3}
    sorted_genes = sorted(v2_dict.items(), key=lambda x: (
        category_order.get(x[1].get("category", "unknown"), 99),
        x[0],
    ))

    for gene_key, info in sorted_genes:
        output["gene_types"][gene_key] = {
            "category": info["category"],
            "scope": info["scope"],
            "synonyms": info["synonyms"],
        }

    with open(output_path, 'w', encoding='utf-8') as f:
        if HAS_YAML:
            yaml.dump(output, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        else:
            # Manual YAML output as fallback
            f.write(f"version: \"2.0\"\n")
            f.write(f"metadata:\n")
            f.write(f"  source: \"g2t v2.0 + SynGenes integration\"\n")
            f.write(f"  description: \"Unified gene dictionary with organism-aware scope\"\n")
            f.write(f"gene_types:\n")
            for gene_key, info in output["gene_types"].items():
                f.write(f"  {gene_key}:\n")
                f.write(f"    category: {info['category']}\n")
                f.write(f"    scope: {info['scope']}\n")
                f.write(f"    synonyms:\n")
                for syn in info["synonyms"]:
                    # Quote strings that contain special YAML characters
                    if any(c in syn for c in ':{}[],&*?|->!%@`#'):
                        f.write(f'      - "{syn}"\n')
                    else:
                        f.write(f"      - {syn}\n")

    return output_path


def print_stats(v2_dict):
    """Print statistics about the v2 dictionary."""
    by_category = defaultdict(lambda: {"count": 0, "synonyms": 0})
    for gene_key, info in v2_dict.items():
        cat = info["category"]
        by_category[cat]["count"] += 1
        by_category[cat]["synonyms"] += len(info["synonyms"])

    print("\n=== Gene Dictionary v2.0 Statistics ===\n")
    total_genes = 0
    total_synonyms = 0
    for cat in ["mitochondrial", "chloroplast", "nuclear", "bacterial", "unknown"]:
        if cat in by_category:
            stats = by_category[cat]
            print(f"  {cat:15s}: {stats['count']:4d} genes, {stats['synonyms']:5d} synonyms")
            total_genes += stats["count"]
            total_synonyms += stats["synonyms"]

    print(f"\n  {'TOTAL':15s}: {total_genes:4d} genes, {total_synonyms:5d} synonyms")


def main():
    parser = argparse.ArgumentParser(description="Build g2t v2.0 gene dictionary")
    parser.add_argument("--output", "-o", default=None, help="Output YAML path")
    parser.add_argument("--dry-run", action="store_true", help="Print stats only, don't write file")
    args = parser.parse_args()

    g2t_root = find_g2t_root()
    syngenes_root = find_syngenes_root()

    # Paths
    g2t_yaml = g2t_root / "src" / "g2t" / "data" / "gene_dict.yaml"
    mt_dir = syngenes_root / "data Names" / "mt"
    cp_dir = syngenes_root / "data Names" / "cp"

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = g2t_root / "src" / "g2t" / "data" / "gene_dict_v2.yaml"

    print(f"Reading g2t v1 dictionary: {g2t_yaml}")
    print(f"Reading SynGenes mt data:  {mt_dir}")
    print(f"Reading SynGenes cp data:  {cp_dir}")

    # Read existing g2t dictionary
    g2t_dict = read_g2t_gene_dict(g2t_yaml)
    print(f"  Found {len(g2t_dict)} existing gene types")

    # Build v2 dictionary
    v2_dict = build_v2_dict(g2t_dict, mt_dir, cp_dir)

    # Print stats
    print_stats(v2_dict)

    if args.dry_run:
        print("\n[DRY RUN] No file written.")
        return

    # Write output
    write_v2_yaml(v2_dict, output_path)
    print(f"\nWritten to: {output_path}")


if __name__ == "__main__":
    main()
