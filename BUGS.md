# g2t 已知问题与解决方案

本文档记录 g2t 使用过程中发现的问题及其解决方案。

---

## 未匹配序列分析 (2026-05-12)

**功能**: g2t 现在会自动生成未匹配序列报告，帮助用户了解哪些序列未被分类。

**报告位置**: `labeled_genes/unmatched_report.txt`

**常见未匹配类型**:

| 类别 | 说明 | 处理建议 |
|------|------|----------|
| mRNA/cDNA | cDNA 克隆序列 | 正常忽略，非标准基因标记 |
| TPA_assembly | 第三方注释组装数据 | 正常忽略 |
| whole_genome_shotgun | 基因组鸟枪法测序 | 正常忽略 |
| genome_assembly | 基因组组装数据 | 正常忽略 |
| creatine_kinase | 肌酸激酶基因 | 可添加到 gene_dict.yaml |
| ATP_synthase | ATP 合酶基因 | 可添加到 gene_dict.yaml |

**如何添加新基因类型**:
编辑 `src/g2t/data/gene_dict.yaml`，添加新的基因类型和同义词。

---

## 2026-05-12: Exit code 1 但处理成功

**问题描述**:
运行 `g2t` 命令时，即使所有步骤都成功完成，程序仍返回 exit code 1。

**原因**:
当存在未匹配序列时（如基因组数据、TPA_asm 数据），程序会输出警告但继续处理。某些情况下这会导致非零退出码。

**解决方案**:
检查输出文件是否完整生成。如果 `organized_species_voucher.csv` 存在且有内容，则处理成功。

**验证方法**:
```bash
# 检查最终输出文件
ls -la output/organized_genes/organized_species_voucher.csv
```

---

## 2026-05-12: TaxonID 查询失败

**问题描述**:
日志中出现 `WARNING - TaxonID '3268408' query error: 3268408 taxid not found`

**原因**:
ete3 的 NCBI taxonomy 数据库可能未更新，或该 TaxonID 是新提交的序列尚未收录。

**影响**:
仅影响 Class/Order/Family/Genus 列的自动填充，不影响基因分类和凭证构建。

**解决方案**:
- 忽略此警告（不影响核心功能）
- 或更新 ete3 数据库: `ete3 -u`

---

## 2026-05-12: 未匹配序列数量较多

**问题描述**:
处理后有较多未匹配序列（如 426 个）。

**原因**:
这些通常是：
1. 基因组组装数据 (whole genome shotgun)
2. TPA_asm 数据 (Third Party Annotation)
3. 非标准基因命名

**解决方案**:
这些序列不是标准的基因序列，正常情况下可以忽略。如需处理：
1. 检查 `unmatched_sequences.csv` 了解具体内容
2. 在 `gene_dict.yaml` 中添加新的同义词

---

## 性能优化记录

| 日期 | 优化项 | 效果 |
|------|--------|------|
| 2026-05-12 | organize.py: iterrows() → itertuples() | 50-100x 加速 |
| 2026-05-12 | classify.py: to_dict() → itertuples() | 14% 加速 |
| 2026-05-12 | extract.py: 移除双重文件读取 | 19% 加速 |

---

## 提交记录

| Commit | 日期 | 描述 |
|--------|------|------|
| c73bd81 | 2026-05-12 | 性能优化 + cob 同义词修复 |
| 6c6c186 | 2026-05-12 | Initial release v1.0.0 |

---

## 已知限制（v0.01，有意保留）：跨区段 rRNA 记录按整条归类

**现象**：classify 按 DEFINITION 把整条记录归入一个基因类型。同时包含多个区段的 rRNA 记录被归入 `its1-its2`，其中的 18S 或 28S 序列不会出现在 18s/28s 列。

**实例（Priapulidae，175 条 rRNA 记录中 9 条）**：
- AY210840（*Priapulus caudatus*，4271 bp）、AH010828（*Halicryptus spinulosus*，4135 bp）：5.8S + ITS2 + 约 3.66 kb 28S，特征表有区段坐标。
- OP247717、OP247728、PQ326428 等 7 条：18S（部分）+ ITS1，1.4–1.8 kb，特征表只有一个笼统的 misc_RNA，无区段边界。

**决定**：暂按原逻辑。备选方案：区段层拆分——特征表有坐标的直接切；无坐标的用 ITSx/barrnap 识别边界；都不行标"未拆分"。需要时再实施。

**手动处理**：在 Records 或 assigned_genes_types_all.csv 中筛选 gene_type = its1-its2 且长度 > 1000 bp 的记录，核对 DEFINITION。


---

## 已知未修复（v0.02 代码审查遗留，2026-10-09）

审查报告：yzz 工作流 `docs/g2t_review_legacy.md`、`docs/g2t_review_new_code.md`。

1. **第 2 轮复查**（classify `recheck_match_row`）：`Topology=circular` 直接判 mtgenome，不看长度和细胞器；第 2 轮不使用 `--mtgenes_list/--ntgenes_list`；`--which_gene_types_extract` 只含第 1 轮结果（L-18）。
2. **`g2t` 的 `--extract_extra/--classify_extra/--voucher_extra/--organize_extra`**：解析后未传给 `run()`，目前无效；需要改单步参数时请分步运行各命令。
3. **格式损坏的记录**：Biopython 解析失败后，同一文件中其后的记录不再读取；`extraction_report.json` 标为部分完成并给出应有/实际记录数，但不恢复。
4. **矩阵元数据逐列取首个非空值**（organize `first_nonempty`）：地点、日期等可能来自同一标本的不同记录（L-09 只修了 Conflict 列）。
5. **mtgenome / 18-28s 的列复制**不检查记录是否注释了该基因（默认行为保留；可用 `--mtgenome_includes` 限定）。
