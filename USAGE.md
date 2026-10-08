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
# 运行完整流程 (推荐使用 --stream 处理大文件)
g2t -i /path/to/genbank_files -o /path/to/output --stream

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
│   └── assigned_genes_types_all.csv # Step 2 输出：带基因类型标签
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
