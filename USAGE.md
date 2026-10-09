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
- ete3: NCBI 分类学查询（用于 Class/Order/Family/Genus；不可用时——如 Python 3.13——改为按 TaxonID 从 NCBI Taxonomy 在线获取一次并缓存到 `gb_metadata/taxonomy_cache.json`；无网络时这几列为空并给出警告）
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

- 按 accession 列表分批下载（不依赖会过期的 WebEnv），失败自动重试。重跑同一命令即续传：每次都从 NCBI 重新取 accession 列表；批文件只有在恰好包含该批的 accession 时才复用（新增或撤回的记录会被补上或去掉）；上次运行多出的批文件会删除。
- 每种选择存独立子目录，按选项命名（markers / mito / mitogenome / gene-COI_18S / incl-wgs / all-mito / `--query` 时为 custom-<检索式哈希> …，或 `--tag 名称`），附 `accessions.tsv` 与 `manifest.json`（检索式、日期、accession 列表校验值、是否完整）。为另一个检索式建的目录不会被复用（报错，请换 `--tag`）。`--dry-run` 不写任何文件。
- 邮箱与 API key：`-e/-k` 或环境变量 `NCBI_EMAIL` / `NCBI_API_KEY`（可选；有 key 时 10 次/秒）。
- 旧脚本 `scripts/download_entrez.py` 现为 `g2t-download` 的包装，参数兼容，但默认选择（markers，而不是全部记录）、输出位置（`<out>/<tag>/`）和批大小（200，而不是 500）已变；要旧的默认行为请加 `--all`。

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
| 候选 | 同一物种 + 共享核心编号：去掉基因名、末尾注释如 "(holotype)" 和标点；`2014-1234` 这类年份–序号整体保留；`WS399, WS400` 这类列表每个编号各算一个；≥5 位的数字去掉前缀后也可配对（`USNM 123456` ~ `123456`）。不同物种同号、或物种名缺失，只报告、不合并 |
| 强 | 同一篇论文（不含 Direct Submission）、采集日期完全相同、坐标相差 ≤0.01° |
| 中 | 第一作者相同、采集人相同、详细地点相同、坐标相差 ≤0.5° |
| 矛盾（阻止合并） | 日期不符、国家不同、坐标相差 >0.5° |

判定：有强证据 → 合并（high）；中等证据 ≥2 → 合并（medium）；否则不合并。两组已含同一基因时只接受强证据。合并两个标本前会检查它们的每一对记录，所以 A–B–C 这样的链条不会把互相矛盾、或会让同一基因在一个标本中出现两次而又没有强证据的 A 和 C 并在一起。

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

### 矩阵中的凭证号

矩阵在 `organism` 之后有三列：`voucher_as_submitted`（该标本在 GenBank 中出现过的全部 specimen_voucher 原样写法，用于看出提交差异）、`voucher_standardized`（每个标本一个统一凭证号：去基因名前缀、取信息最全的写法，按 INSDC `/specimen_voucher` 格式 `[机构代码:[收藏代码:]]标本号` 书写，如 `COI_ZMMU_MSU_WS2585; 28S_ZMMU_WS2585; WS2585` → `ZMMU:MSU:WS2585`，`ZMMU:WS30980; ZMMU_WS30980` → `ZMMU:WS30980`）、`voucher_note`（几种写法不只是前缀/分隔符不同时提示人工核对）。

### Step 5: 校正元数据 (curate，可选)

#### 目的

序列通常在投稿时提交到 GenBank。审稿后常有改动：标本重新鉴定或组合变更、被描述为新种、题目改变、坐标和日期更正。GenBank 记录只能由提交者更新。INSDC 欢迎"作者对错误的更正和记录的更新"（[INSDC policy](https://www.insdc.org/policy/)）；NCBI 要求文章发表后由提交者把出版信息发给 GenBank（[GenBank 概述](https://www.ncbi.nlm.nih.gov/genbank/)），来源信息按表格提交更新（[Update GenBank records](https://www.ncbi.nlm.nih.gov/genbank/update/)）。实际上很多记录从未更新，下载的数据集里元数据可能过时或错误。公共序列鉴定错误和注释不全已有文献记载（Bridge et al. 2003, *New Phytologist* 160: 43–48, [doi:10.1046/j.1469-8137.2003.00861.x](https://doi.org/10.1046/j.1469-8137.2003.00861.x)；Meiklejohn et al. 2019, *PLOS ONE* 14: e0217084, [doi:10.1371/journal.pone.0217084](https://doi.org/10.1371/journal.pone.0217084)）。

Step 5 让分类学者把核实过的信息用到"标本 × 基因"矩阵上，不修改 GenBank，也不丢失原始值：

- **A**：由 GenBank 整理的矩阵，保持提交时的原样（一个标本一行）。
- **B**：你自己整理的表，沿用你已有的格式。内容是你核实过的信息：有效种名、正式发表的题目和作者、更正后的坐标或日期，以及水深等额外字段。
- **C**：用 B 校正 A 后的新矩阵，用于后续分析。

原则：

- **以标本为单位**：B 的每一行校正一个标本（矩阵中的一行）。B 中的凭证号和登录号只用来找到这一行。
- **保守**：B 的某行匹配不到标本、匹配到多个标本、或凭证号与登录号指向不同标本时，不应用，列出原因。
- **可追溯**：C 中保留 A 的原样。每处修改都记录原值、新值和匹配依据，改动的单元格标黄。修改记录也就是需要告知提交者或 NCBI 的更正清单。
- **沿用 GenBank 字段与格式**：值按 GenBank 格式写（`lat_lon` `66.55 N 33.10 E`，`collection_date` `12-Jun-2019`），使 C 与新下载的数据可直接比较。Darwin Core 字段命名和导出 NCBI 来源信息更新表计划后续加入。

g2t 不修改 GenBank。C 是你的校正副本；INSDC 中的记录仍由提交者负责。

#### 用法

NCBI 矩阵（**A**）+ 你自己整理的表（**B**）→ 校正后的新矩阵（**C**）。

```bash
g2t-curate -m matrix_A.xlsx -u 我的表_B.xlsx -o 校正后_C.xlsx
g2t-curate -m matrix_A.xlsx -u B.xlsx -o C.xlsx --map "编号=voucher" --map "Lat=latitude"   # 识别错时强制指定
```

**B 可以是你已有的任意表格**（Excel/csv，列名随意，中英文均可；`--sheet` 选工作表）。自动识别：

| B 中的列（示例） | 用途 |
|---|---|
| voucher、凭证号、标本号、catalog number | 定位标本（任何写法） |
| 任何登录号列：`COI登录号`、`28S accession`，或值像 `ON792938.1` 的列 | 定位标本；同一行的登录号必须指向同一标本 |
| species、拉丁名、生物种拉丁名 | `organism` |
| 纬度 + 经度（十进制） | 合并为 GenBank 格式 `lat_lon`（`66.55 N 33.10 E`） |
| 采集日期（`20190612`、`2019-06-12`） | `collection_date`，GenBank 格式（`12-Jun-2019`） |
| 采集地、国家、论文题目、作者、期刊、采集人、鉴定人 | `geo_loc_name`、`country`、`Ref1Title`、`Ref1Authors`、`Ref1Journal`、`collected_by`、`identified_by` |
| 其余列（如 `水深(m)`） | 作为新列 `user:<列名>` 加入 C |

C 的工作表：*Matrix*（改过的单元格标黄）、*Changes*（原值 → 新值、匹配依据）、*Problems*（匹配不到、匹配到多个标本、登录号与凭证号指向不同标本——均不应用）、*Column mapping*（B 每列如何使用）、*Matrix (GenBank)*（A 原样）。与 A 相同的值不算修改。

也可用预填模板：

```bash
g2t-curate -m matrix.csv --template corrections.xlsx      # 生成模板：每个标本一行，已填当前值，另存 baseline 工作表
# 在 corrections.xlsx 中改需要校正的单元格；空单元格 = 不改
g2t-curate -m matrix.csv -u corrections.xlsx -o curated.csv
```

- 定位：`accession`（任一基因列的登录号，版本号可省）和/或 `voucher`（任何写法，如 `WS2585` 可找到 `ZMMU_MSU_WS2585`）；同一凭证号出现在多个物种时加 `organism_match`。其余列为要设置的值（`organism`、`lat_lon`、`geo_loc_name`、`collection_date`、`Ref1Title`、`curation_source` 等，也可是新列）。
- 输出：`curated.csv`（含 `curated_fields` 列）、`curated_curation_log.csv`（每处修改：原值 → 新值、匹配依据）、`curated_curation_problems.csv`（匹配不到、匹配到多个标本、accession 与 voucher 指向不同标本的行——均不应用）。原矩阵不改。
- 使用模板时只有你改过的单元格才算修改（与模板中的 `baseline (do not edit)` 工作表比较），没改的单元格不会把 GenBank 之后更新的值改回去；你改过而 GenBank 也已变化的单元格列为冲突，不应用。
- 不能通过校正把单元格清空（空 = 不改）。

---

## 输出文件说明

```
output/
├── gb_metadata/
│   └── final.csv                    # Step 1 输出：所有序列的元数据
├── labeled_genes/
│   ├── assigned_genes_types_all.csv # Step 2 输出：带基因类型标签
│   ├── unmatched_sequences.csv      # Step 2：未识别出基因类型的记录
│   ├── filtered_records.csv         # Step 2：被过滤掉的记录及原因（长度、缺长度、重复 LocusID 等）
│   └── record_status.csv            # Step 2：每条输入记录一行：assigned / unmatched / filtered + 原因
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
- **断点续传**: 使用 `--resume` 参数跳过已完成步骤；凭证号核对的设置与上次不同时（记录在 `OUT/.g2t_params.json`），Step 3b 和 4 会重跑。同一目录重跑 `g2t-classify` 会先清除该步骤上次的输出。

---

## 不足与已知限制

1. **基因类型按 DEFINITION 整条判定。** 一条记录只得到一个类型。跨区段的 rRNA 记录（如"5.8S … ITS2 … 28S""18S … ITS1"）被归为 `its1-its2`，其中的 18S 或 28S 部分不会出现在 18s/28s 列（Priapulidae：175 条 rRNA 记录中 9 条，如 AY210840 含约 3.7 kb 28S）。详见 [BUGS.md](BUGS.md)。
2. **不切分序列。** 线粒体基因组（`mtgenome`）的登录号会被复制到全部线粒体基因列（coi、16s、12s、cob、cox2、cox3），18S–ITS–28S 记录（`18-28s`）复制到 18s、28s、its1-its2，**不论该记录是否真的注释了这个基因**（如 FN689349 没有 12S 注释，仍出现在 12s 列）。可用 `g2t-organize --mtgenome_includes coi 16s` 限定列。不会切出各基因的序列；比对前需另行提取。
3. **长度过滤。** 默认长度 150–50,000 bp 以外（以及缺长度、LocusID 重复）的记录在分类前被去掉。它们不会静默消失：逐条列在 `filtered_records.csv`（含原因），`record_status.csv` 给每条输入记录恰好一个状态。若要让它们也参与分类，三个范围都要放宽：`g2t-classify … --length_range2_all 1:1000000 --length2_mtgenes 1:1000000 --length2_ntgenes 1:1000000`。
4. **只识别 13 类基因，按关键词匹配。** 其他位点（核蛋白编码基因、微卫星、Hox 基因等）进入 `unmatched_sequences.csv`，除非在 `gene_dict.yaml` 中添加。DEFINITION 中的标本编号（`isolate CO2`、`voucher COI-12`）和英文单词 "its" 不参与匹配，并排除了几种形似的情况（16S rRNA methyltransferase、histone H3 lysine/demethylase、elongation factor-1 beta/gamma），但写法特殊时关键词匹配仍可能出错。
5. **凭证号与元数据。** Step 3 把同一物种中凭证号字符串完全相同的记录直接合并，不做进一步检查；不同标本恰好同号时会被误并。Step 3b 依赖提交到 GenBank 的元数据：没有共同的论文、日期、坐标或采集人时，同一凭证号的不同写法不会合并（宁缺毋滥）。矩阵中的地点、日期等元数据逐列取该标本第一条有值的记录，因此可能来自不同记录；任一记录被标记冲突，`Conflict` 即为 True。
6. **物种名按提交原样。** 不与 WoRMS、NCBI 异名等分类权威核对；错误鉴定或过时名称保持原样，除非在 Step 5 中校正。
7. **下载。** 只查 NCBI nuccore（不含 BOLD、仅在 ENA 的数据、SRA）。`--gene` 按基因字段和标题关键词检索，写法特殊的记录可能漏检；按基因检索也会返回含该基因的线粒体基因组。
8. **格式损坏的文件。** GenBank 文件中某条记录无法被 Biopython 解析时，该文件中其后的记录不会被读取；`extraction_report.json` 中该文件标为部分完成，并在警告中给出应有与实际解析的记录数。
9. **校正（Step 5）。** B 只校正 A 中已有的标本。没有 GenBank 记录的标本列为"匹配不到"，不会加入 C。空单元格表示"不改"，因此不能通过 B 清空 A 中的值。列识别是根据列名和"像登录号的值"的启发式判断，请查看 *Column mapping* 表，识别错时用 `--map` 指定。B 中的值被视为正确：程序不会把 B 本身与文献或 WoRMS 核对。
10. **代码状态。** v0.02 为早期版本，在 Python 3.13 上以 344 项单元测试验证；早期模块仍有 `ruff` 代码风格警告。

---

## 引用

g2t 目前没有发表论文，请引用软件本身及所用版本：

> Yang, D. (2026). *g2t: GenBank to Taxonomy* (version v0.02) [Computer software]. GitHub. https://github.com/deyuanyang92-dev/gb2taxonomy

```bibtex
@software{yang_g2t_2026,
  author  = {Yang, Deyuan},
  title   = {g2t: GenBank to Taxonomy},
  version = {v0.02},
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
