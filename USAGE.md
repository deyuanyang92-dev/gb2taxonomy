# g2t 使用说明

## 软件简介

g2t (GenBank to Taxonomy) 是一个生物信息学工具，用于将来自同一标本的不同分子数据整理成矩阵格式，便于后续的系统发育分析。

在分类学和系统发育学研究中，我们通常需要对同一标本测序多个基因标记（如 COI、16S、18S 等）。g2t 可以自动从 GenBank 文件中提取这些信息，识别基因类型，通过标本凭证号关联序列，最终生成物种 × 基因的矩阵表。

## 依赖要求

**核心依赖** (必需):
- Python ≥ 3.9
- pandas ≥ 1.5
- biopython ≥ 1.80

**可选依赖**:
- pyyaml: 支持自定义基因词典
- ete3: NCBI 分类学查询
- pytest: 开发测试

## 安装

### 方法 1: 使用现有环境 (推荐)

如果你的环境已经有 pandas 和 biopython，可以直接安装：

```bash
pip install -e /path/to/gb2taxonomy
```

### 方法 2: 创建新环境

```bash
# 使用 mamba (推荐)
mamba create -n g2t python=3.10 -y
mamba activate g2t
pip install -e /path/to/gb2taxonomy

# 或使用 conda
conda create -n g2t python=3.10 -y
conda activate g2t
pip install -e /path/to/gb2taxonomy

# 或使用 pip + venv
python -m venv g2t_env
source g2t_env/bin/activate  # Linux/Mac
# g2t_env\Scripts\activate   # Windows
pip install -e /path/to/gb2taxonomy
```

### 方法 3: 仅安装依赖

```bash
pip install pandas biopython
```

### 验证安装

```bash
g2t --help
```

如果缺少依赖，程序会自动提示安装方法。

---

## 快速开始

```bash
# 0. 先看 NCBI 中该类群有哪些记录, 再下载标记基因记录
g2t-download -t Priapulidae -o gb --dry-run
g2t-download -t Priapulidae -o gb                      # -> gb/markers/batch_0001.gb ...

# 1–4. 运行完整流程 (推荐使用 --stream 处理大文件)
g2t -i gb/markers -o out --stream
# 最终矩阵: out/organized_genes/organized_species_voucher.csv

# 使用 --resume 跳过已完成的步骤
g2t -i /path/to/files -o /path/to/output --resume
```

---

## 处理流程

### Step 0: 下载 (download)

按类群从 NCBI 下载 GenBank 记录。默认只下"标记基因"：排除 WGS contig、mRNA、RefSeq（预测模型及 NC_ 等与 INSDC 重复的拷贝），大类群中这几类常占 95% 以上，对多基因表无用。

```bash
g2t-download -t Priapulidae -o gb_out --dry-run                # 类群构成 + 各基因条数 + 检索式
g2t-download -t Priapulidae -o gb_out                          # gb_out/markers/
g2t-download -t Priapulidae -o gb_out --mito                   # 线粒体记录
g2t-download -t Priapulidae -o gb_out --mitogenome             # 完整线粒体基因组 (10–30 kb)
g2t-download -t Priapulidae -o gb_out --gene COI,18S,28S       # 指定基因 (g2t/ncbi_genes.py)
g2t-download -t Priapulidae -o gb_out --gene COI --minlen 500
g2t-download -t Priapulidae -o gb_out --query 'Russia[Country]'
g2t-download -t Priapulidae -o gb_out --all                    # 全部 (可能极大)
```

- 按 accession 列表分批下载（不依赖会过期的 WebEnv），每批校验条数，失败自动重试；中断后重跑同一命令即续传。
- 每种选择存独立子目录（markers / mito / mitogenome / gene-COI_18S …），附 `accessions.tsv` 与 `manifest.json`（检索式、日期、条数）。
- 邮箱与 API key：`-e/-k` 或环境变量 `NCBI_EMAIL` / `NCBI_API_KEY`（可选；有 key 时 10 次/秒）。
- 旧脚本 `scripts/download_entrez.py` 仍可用，参数相同。

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

### Step 3b: 凭证号核对 (reconcile)

同一标本的不同基因，提交者常把凭证号写成不同形式（`COI_ZMMU_MSU_WS399` / `28S_ZMMU_WS399` / `WS399`），Step 3 按字面匹配会拆成多行。本步骤依据记录中的证据合并：

| 证据 | 内容 |
|---|---|
| 候选 | 同一物种 + 核心编号相同（去掉基因名与馆藏缩写，≥4 字符）；不同物种同号只报告、不合并 |
| 强 | 同一篇论文（不含 Direct Submission）、采集日期完全相同、坐标相差 ≤0.01° |
| 中 | 第一作者相同、采集人相同、详细地点相同、坐标相差 ≤0.5° |
| 矛盾（阻止合并） | 日期不符、国家不同、坐标相差 >0.5° |

判定：有强证据 → 合并（high）；中等证据 ≥2 → 合并（medium）；否则不合并。两组已含同一基因时只接受强证据。

```bash
g2t-reconcile -i updated_species_voucher.csv -o /path/to/out [--min_confidence high]
```

输出 `reconciled_species_voucher.csv`（新增 species_voucher_g2t / voucher_core / match_basis / match_confidence / match_evidence）与 `reconcile_report.csv`（每对候选的证据与判定）。完整流程中默认开启；`g2t --skip_reconcile` 关闭，`--reconcile_min_confidence high` 只按强证据合并。

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
│   ├── assigned_genes_types_all.csv # Step 2 输出：带基因类型标签
│   └── unmatched_sequences.csv      # Step 2：未识别出基因类型的记录
├── updated_species_vouchers/
│   ├── updated_species_voucher.csv  # Step 3 输出：带标本凭证号
│   ├── reconciled_species_voucher.csv # Step 3b 输出：凭证号核对后
│   └── reconcile_report.csv         # Step 3b：每对候选的证据与判定
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

### Q: 提示缺少依赖怎么办？

程序会自动检测并提示安装方法。如果看到错误信息：

```
ERROR: Missing required dependencies!
Missing: pandas, biopython
```

请按照提示安装：

```bash
pip install pandas biopython
```

---

## 性能建议

- **大文件处理**: 使用 `--stream` 参数，避免内存溢出
- **多文件处理**: 使用 `--batch` 参数并行处理
- **断点续传**: 使用 `--resume` 参数跳过已完成步骤

---

## 不足与已知限制

1. **基因类型按 DEFINITION 整条判定。** 一条记录只得到一个类型。跨区段的 rRNA 记录（如"5.8S … ITS2 … 28S""18S … ITS1"）被归为 `its1-its2`，其中的 18S 或 28S 部分不会出现在 18s/28s 列（Priapulidae：175 条 rRNA 记录中 9 条，如 AY210840 含约 3.7 kb 28S）。详见 [BUGS.md](BUGS.md)。
2. **不切分序列。** 线粒体基因组（`mtgenome`）或 18S–ITS–28S 记录（`18-28s`）的登录号会被复制到它覆盖的各基因列，但不会切出各基因的序列；比对前需另行提取。
3. **长度过滤。** 默认长度 150–50,000 bp 以外的记录在分类前被去掉，且**不会**出现在 `unmatched_sequences.csv`。要列出它们：`g2t-classify … --length_range2_all 1:1000000`。
4. **只识别 13 类基因。** 其他位点（核蛋白编码基因、微卫星、Hox 基因等）进入 `unmatched_sequences.csv`，除非在 `gene_dict.yaml` 中添加。
5. **凭证号。** Step 3 把同一物种中凭证号字符串完全相同的记录直接合并，不做进一步检查；不同标本恰好同号时会被误并。Step 3b 依赖提交到 GenBank 的元数据：没有共同的论文、日期、坐标或采集人时，同一凭证号的不同写法不会合并（宁缺毋滥）。
6. **物种名按提交原样。** 不与 WoRMS、NCBI 异名等分类权威核对；错误鉴定或过时名称保持原样。
7. **下载。** 只查 NCBI nuccore（不含 BOLD、仅在 ENA 的数据、SRA）。`--gene` 按基因字段和标题关键词检索，写法特殊的记录可能漏检；按基因检索也会返回含该基因的线粒体基因组。
8. **代码状态。** v0.01 为早期版本，在 Python 3.13 上以 232 项单元测试验证；早期模块仍有 `ruff` 代码风格警告。

---

## 引用

g2t 目前没有发表论文，请引用软件本身及所用版本：

> Yang, D. (2026). *g2t: GenBank to Taxonomy* (version v0.01) [Computer software]. GitHub. https://github.com/deyuanyang92-dev/gb2taxonomy

```bibtex
@software{yang_g2t_2026,
  author  = {Yang, Deyuan},
  title   = {g2t: GenBank to Taxonomy},
  version = {v0.01},
  year    = {2026},
  url     = {https://github.com/deyuanyang92-dev/gb2taxonomy}
}
```

GitHub 页面右侧的 "Cite this repository" 读取 [CITATION.cff](CITATION.cff)。另请引用：

- **序列数据**：所用 GenBank 记录的原始文献（输出中的 `Ref1Authors`、`Ref1Title`、`Ref1Journal` 列）以及 NCBI GenBank。
- **Biopython**（用于下载和解析记录）：Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

---

## 许可证

MIT License
