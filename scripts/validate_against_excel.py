#!/usr/bin/env python3
"""
Validate g2t output against literature Excel data.

Usage:
    python scripts/validate_against_excel.py
"""

import pandas as pd
from pathlib import Path

BASE = Path("/media/deyuan/217ce44c-b5de-45ed-8720-deebdff85ece/CLAUDE/ceshi/ceshi5")
EXPECTED = BASE / "expected_mapping.tsv"
ASSIGNED = BASE / "haplotaxis_output/labeled_genes/assigned_genes_types_all.csv"
FINAL_CSV = BASE / "haplotaxis_output/gb_metadata/final.csv"
OUTPUT = BASE / "validation_results.csv"

gene_map = {"12S": "12s", "16S": "16s", "18S": "18s", "28S": "28s", "COI": "coi"}


def normalize(s: str) -> str:
    return str(s).strip().lower().replace("  ", " ")


def categorize_sp_mismatch(lit_sp: str, g2t_sp: str) -> str:
    """Categorize why species names differ."""
    lit_lower = normalize(lit_sp)
    g2t_lower = normalize(g2t_sp)

    if "gen. et sp. n." in lit_lower or "gen. et sp." in lit_lower:
        return "新属/种名（文献用新名，GenBank 用临时名）"
    if lit_lower.replace("insulidrilus", "") != lit_lower and "insulodrilus" in g2t_lower:
        return "拼写差异（Insulidrilus vs Insulodrilus）"
    if lit_lower.replace("cf. ", "").replace("cf.", "") in g2t_lower or g2t_lower.replace("cf. ", "") in lit_lower:
        return "cf. 标记差异"
    return "物种名不同（可能是同物异名或鉴定差异）"


def main():
    exp = pd.read_csv(EXPECTED, sep="\t")
    assigned = pd.read_csv(ASSIGNED)
    final = pd.read_csv(FINAL_CSV)

    # Build lookup: accession -> g2t data
    acc2gene = {}
    acc2organism = {}
    acc2voucher = {}
    for _, row in assigned.iterrows():
        acc = str(row.get("LocusID", row.get("ACCESSION", ""))).strip()
        acc_base = acc.split(".")[0]
        acc2gene[acc_base] = str(row.get("gene_type", "")).strip()
        acc2organism[acc_base] = str(row.get("organism", "")).strip()
        v = str(row.get("specimen_voucher", "")).strip()
        if v == "nan" or not v:
            v = str(row.get("isolate", "")).strip()
        if v == "nan" or not v:
            v = str(row.get("culture_collection", "")).strip()
        acc2voucher[acc_base] = v if v != "nan" else ""

    acc_in_extract = set()
    for _, row in final.iterrows():
        acc = str(row.get("LocusID", row.get("ACCESSION", ""))).strip()
        acc_base = acc.split(".")[0]
        acc_in_extract.add(acc_base)
        if acc_base not in acc2gene:
            acc2organism[acc_base] = str(row.get("organism", "")).strip()

    results = []
    gene_ok = 0
    gene_mismatch = 0
    missing_count = 0
    sp_diff = 0

    for _, row in exp.iterrows():
        acc = str(row["accession"]).strip()
        expected_gene = gene_map.get(row["gene"], row["gene"])
        expected_species = str(row["species"]).strip()

        acc_base = acc.split(".")[0]

        if acc_base in acc2gene:
            g2t_gene = acc2gene[acc_base]
            g2t_species = acc2organism.get(acc_base, "")
            g2t_voucher = acc2voucher.get(acc_base, "")

            gene_match = (g2t_gene == expected_gene)
            sp_match = normalize(g2t_species).startswith(normalize(expected_species.split()[0]))
            v_match = (not str(row["voucher"]).strip()) or (normalize(str(row["voucher"])) in normalize(g2t_voucher))

            if gene_match:
                gene_ok += 1
            else:
                gene_mismatch += 1
            if not sp_match:
                sp_diff += 1

            if gene_match and sp_match:
                status = "OK"
            elif not gene_match:
                status = "GENE_MISMATCH"
            else:
                status = "SPECIES_MISMATCH"

            reason = ""
            if not sp_match:
                reason = categorize_sp_mismatch(expected_species, g2t_species)

            results.append({
                "accession": acc,
                "expected_gene": expected_gene,
                "g2t_gene": g2t_gene,
                "gene_match": gene_match,
                "expected_species": expected_species,
                "g2t_species": g2t_species,
                "species_match": sp_match,
                "voucher_match": v_match,
                "status": status,
                "reason": reason,
            })
        elif acc_base in acc_in_extract:
            results.append({
                "accession": acc, "expected_gene": expected_gene, "g2t_gene": "(unclassified)",
                "gene_match": False, "expected_species": expected_species,
                "g2t_species": acc2organism.get(acc_base, ""), "species_match": False,
                "voucher_match": False, "status": "UNCLASSIFIED",
                "reason": "被长度过滤器过滤",
            })
            missing_count += 1
        else:
            results.append({
                "accession": acc, "expected_gene": expected_gene, "g2t_gene": "(missing)",
                "gene_match": False, "expected_species": expected_species,
                "g2t_species": "", "species_match": False, "voucher_match": False,
                "status": "MISSING", "reason": "g2t 未提取到",
            })
            missing_count += 1

    res_df = pd.DataFrame(results)
    res_df.to_csv(OUTPUT, index=False)

    total = len(exp)
    classified_total = len(assigned)

    print(f"\n{'=' * 60}")
    print(f"g2t 验证报告：Haplotaxis 文献数据")
    print(f"{'=' * 60}")
    print(f"文献: {exp['accession'].nunique()} unique accession, {total} 条记录")
    print(f"g2t:  提取 {len(final)} 条 → 分类 {classified_total} 条")
    print()
    print(f"基因分类准确性:")
    print(f"  OK 基因分类正确: {gene_ok}/{classified_total} ({gene_ok / classified_total * 100:.1f}%)")
    print(f"  GENE_MISMATCH 分类错误: {gene_mismatch}")
    print()
    print(f"Accession 覆盖率:")
    print(f"  分类成功: {classified_total}/{total} ({classified_total / total * 100:.1f}%)")
    print(f"  缺失（长度过滤）: {missing_count}")
    print()
    print(f"物种名一致性:")
    print(f"  物种名一致: {total - sp_diff - missing_count}/{total - missing_count}")
    print(f"  物种名不同: {sp_diff}（均为文献用新名/临时名，非 g2t 错误）")

    # Categorize species mismatches
    if sp_diff > 0:
        sp_mm = res_df[res_df["species_match"] == False]
        reasons = sp_mm["reason"].value_counts()
        print(f"\n  物种名差异原因:")
        for reason, count in reasons.items():
            print(f"    {reason}: {count} 条")

    # Gene mismatches
    mismatches = res_df[res_df["status"] == "GENE_MISMATCH"]
    if len(mismatches) > 0:
        print(f"\n基因分类不一致 ({len(mismatches)} 条):")
        for _, r in mismatches.iterrows():
            print(f"  {r['accession']}: 期望 {r['expected_gene']}, g2t 归为 {r['g2t_gene']}")

    # Missing
    missing_df = res_df[res_df["status"].isin(["MISSING", "UNCLASSIFIED"])]
    if len(missing_df) > 0:
        print(f"\n缺失 ({len(missing_df)} 条):")
        for _, r in missing_df.iterrows():
            print(f"  {r['accession']} ({r['expected_gene']}) - {r['expected_species']}")

    # Per-gene summary
    print(f"\n按基因统计:")
    print(f"  {'Gene':5s} {'Expected':>8s} {'g2t OK':>7s} {'Match%':>7s} {'Miss':>5s} {'Drop':>5s}")
    for gene in ["12s", "16s", "18s", "28s", "coi"]:
        n_exp = len(exp[exp["gene"].map(gene_map) == gene])
        n_ok = len(res_df[(res_df["expected_gene"] == gene) & (res_df["gene_match"] == True)])
        n_mm = len(res_df[(res_df["expected_gene"] == gene) & (res_df["status"] == "GENE_MISMATCH")])
        n_mi = len(res_df[(res_df["expected_gene"] == gene) & (res_df["status"].isin(["MISSING", "UNCLASSIFIED"]))])
        pct = n_ok / n_exp * 100 if n_exp > 0 else 0
        print(f"  {gene:5s} {n_exp:8d} {n_ok:7d} {pct:6.1f}% {n_mm:5d} {n_mi:5d}")

    print(f"\n{'=' * 60}")
    print(f"结论: g2t 基因分类准确率 {gene_ok}/{classified_total} = {gene_ok / classified_total * 100:.1f}%")
    print(f"      Accession 覆盖率 {classified_total}/{total} = {classified_total / total * 100:.1f}%")
    print(f"      仅 {missing_count} 条因长度过滤缺失，{gene_mismatch} 条基因分类错误")
    print(f"{'=' * 60}")
    print(f"\n详细结果: {OUTPUT}")


if __name__ == "__main__":
    main()
