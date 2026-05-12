# g2t 使用说明

## 软件简介

g2t (GenBank to Taxonomy) 是一个生物信息学工具，用于将来自同一标本的不同分子数据整理成矩阵格式，便于后续的系统发育分析。

在分类学和系统发育学研究中，我们通常需要对同一标本测序多个基因标记（如 COI、16S、18S 等）。g2t 可以自动从 GenBank 文件中提取这些信息，识别基因类型，通过标本凭证号关联序列，最终生成物种 × 基因的矩阵表。

## 安装

```bash
pip install -e .
pip install -e ".[dev]"    # 包含测试依赖
pip install -e ".[yaml]"    # 支持自定义基因词典
```

**依赖**: Python ≥3.9, pandas, biopython

---

## 快速开始

```bash
# 运行完整流程
g2t -i /path/to/genbank_files -o /path/to/output --stream

# 使用 --resume 跳过已完成的步骤
g2t -i /path/to/files -o /path/to/output --resume
```

---

## 处理流程

### Step 1: 提取元数据 (extract)

从 GenBank 文件中提取序列信息和标本信息。

```bash
g2t-extract -i /path/to/gb_files -o /path/to/out --stream
```

**提取内容**:
- 序列信息: LocusID, ACCESSION, 长度, Definition
- 标本信息: organism, specimen_voucher, isolate, strain, country, collected_by
- 分类信息: TaxonID, Class, Order, Family, Genus
- 组装信息: Assembly Method, Sequencing Technology

**参数说明**:
- `--stream`: 流式处理，适合大文件（推荐）
- `--batch`: 并行处理多个文件

**输出文件**: `final.csv`

---

### Step 2: 基因类型识别 (classify)

根据序列的 Definition 字段自动识别基因类型。

```bash
g2t-classify -i final.csv -o /path/to/out
```

**支持的基因类型**:

| 类型 | 基因标记 |
|------|----------|
| 线粒体基因 | coi, 16s, 12s, cob, mtgenome, cox3, cox2 |
| 核基因 | 18-28s, its1-its2, ef-1, 28s, 18s, h3 |

**输出文件**: `assigned_genes_types_all.csv`

---

### Step 3: 构建标本凭证 (voucher)

通过标本凭证号将来自同一标本的不同序列关联起来。

```bash
g2t-voucher -i assigned_genes_types_all.csv -o /path/to/out
```

**凭证号优先级**:
1. specimen_voucher
2. isolate
3. culture_collection
4. clone
5. strain
6. ACCESSION（兜底）

**输出文件**: `updated_species_voucher.csv`

---

### Step 4: 生成矩阵 (organize)

生成物种 × 基因的矩阵表，每行一个标本，每列一个基因标记。

```bash
g2t-organize -i updated_species_voucher.csv -o output.csv
```

**输出格式**:

| species_voucher_new | organism | coi | 16s | 18s | ... |
|---------------------|----------|-----|-----|-----|-----|
| SpeciesA_voucher1 | Species A | LC123456 | LC123457 | LC123458 | |
| SpeciesB_voucher2 | Species B | LC123459 | | LC123460 | |

单元格内为该标本对应基因的 LocusID，多个序列用分号分隔。

**输出文件**: `organized_species_voucher.csv`

---

## 输出文件说明

```
output/
├── gb_metadata/
│   └── final.csv                    # Step 1 输出：所有序列的元数据
├── labeled_genes/
│   └── assigned_genes_types_all.csv # Step 2 输出：带基因类型标签
├── updated_species_vouchers/
│   └── updated_species_voucher.csv  # Step 3 输出：带标本凭证号
└── organized_genes/
    └── organized_species_voucher.csv # Step 4 输出：最终矩阵
```

---

## 常见问题

### Q: 如何处理没有 specimen_voucher 的序列？

程序会依次尝试 isolate → culture_collection → clone → strain，最终使用 ACCESSION 作为兜底。

### Q: 如何添加自定义基因类型？

编辑 `src/g2t/data/gene_dict.yaml`，添加新的基因类型和同义词：

```yaml
gene_types:
  new_gene:
    synonyms:
      - synonym 1
      - synonym 2
```

或使用外部词典：

```bash
g2t-classify -i input.csv -o output --path_dict custom_dict.yaml
```

### Q: 如何只运行部分步骤？

```bash
g2t -i input -o output --skip_extract    # 跳过 Step 1
g2t -i input -o output --skip_classify   # 跳过 Step 2
```

---

## 引用

如果在研究中使用 g2t，请引用：

```bibtex
@article{g2t2026,
  title = {g2t: a Python pipeline for GenBank-to-Taxonomy gene type classification and species voucher organization},
  journal = {Bioinformatics},
  year = {2026},
}
```

---

## 许可证

MIT License
