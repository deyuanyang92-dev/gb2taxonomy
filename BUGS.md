# g2t 已知问题与解决方案

本文档记录 g2t 使用过程中发现的问题及其解决方案。

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
