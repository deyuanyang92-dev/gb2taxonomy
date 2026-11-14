# g2t 使用说明

## 软件简介

g2t (GenBank to Taxonomy) 是一个生物信息学工具，用于将来自同一标本的不同分子数据整理成矩阵格式，便于后续的系统发育分析。

在分类学和系统发育学研究中，我们通常需要对同一标本测序多个基因标记（如 COI、16S、18S 等）。g2t 可以自动从 GenBank 文件中提取这些信息，识别基因类型，通过标本凭证号关联序列，最终生成物种 × 基因的矩阵表。

## 依赖要求

**核心依赖** (必需):
- Python ≥ 3.9
- pandas ≥ 1.5
- biopython ≥ 1.80
- openpyxl ≥ 3.0: Step 6（`g2t-curate`）读写 Excel 文件
- pyyaml ≥ 6.0: 读取自带的 `gene_dict.yaml`，这是同义词的唯一来源

**可选依赖**:
- ete3: NCBI 分类学查询（用于 Class/Order/Family/Genus；不可用时——如 Python 3.13——改为按 TaxonID 从 NCBI Taxonomy 在线获取一次并缓存到 `gb_metadata/taxonomy_cache.json`；无网络时这几列为空并给出警告）
- pytest: 开发测试

## 安装

### 方法 1: 使用现有环境 (推荐)

如果你的环境已经有 pandas、biopython、openpyxl 和 pyyaml，可以直接安装（缺的会由 pip 自动补装）：

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
pip install pandas biopython openpyxl pyyaml
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

图解（GenBank 记录的组成、整理流程图、三张逻辑图）：[docs/genbank-and-workflow.md](docs/genbank-and-workflow.md)。

### Step 0: 下载 (download)

按类群从 NCBI 下载 GenBank 记录。默认只下"标记基因"：排除 WGS contig、mRNA、TSA（转录组组装）contig、RefSeq（预测模型及 NC_ 等与 INSDC 重复的拷贝），大类群中这几类常占 95% 以上，对多基因表无用。

```bash
g2t-download -t Priapulidae -o gb_out --dry-run                # 类群构成（WGS、mRNA、TSA 等）+ 各基因条数 + 检索式
g2t-download -t Priapulidae -o gb_out                          # gb_out/markers/
g2t-download -t Priapulidae -o gb_out --mito                   # 线粒体记录
g2t-download -t Priapulidae -o gb_out --mitogenome             # 完整线粒体基因组 (≥10 kb, 无上限; 含 "genome assembly, organelle: mitochondrion")
g2t-download -t Priapulidae -o gb_out --gene COI,18S,28S       # 指定基因 (g2t/ncbi_genes.py)
g2t-download -t Priapulidae -o gb_out --gene COI --minlen 500
g2t-download -t Priapulidae -o gb_out --query 'Russia[Country]'
g2t-download -t Priapulidae -o gb_out --include-tsa            # 标记基因 + TSA contig (默认以 NOT "tsa"[Keyword] 排除)
g2t-download -t Priapulidae -o gb_out --all                    # 全部 (可能极大)
```

- 按 accession 列表分批下载（不依赖会过期的 WebEnv），失败自动重试。重跑同一命令即续传：每次都从 NCBI 重新取 accession 列表；批文件只有在恰好包含该批的 accession 时才复用（新增或撤回的记录会被补上或去掉）；上次运行多出的批文件会删除。
- 每种选择存独立子目录，按选项命名（markers / mito / mitogenome / gene-COI_18S / incl-wgs / all-mito / `--query` 时为 custom-<检索式哈希> …，或 `--tag 名称`），附 `accessions.tsv` 与 `manifest.json`（检索式、日期、accession 列表校验值、是否完整）。为另一个检索式建的目录不会被复用（报错，请换 `--tag`）。`--dry-run` 不写任何文件。
- 邮箱与 API key：`-e/-k` 或环境变量 `NCBI_EMAIL` / `NCBI_API_KEY`（可选；有 key 时 10 次/秒）。
- 旧脚本 `scripts/download_entrez.py` 现为 `g2t-download` 的包装，参数兼容，但默认选择（markers，而不是全部记录）、输出位置（`<out>/<tag>/`）和批大小（500）已变；要旧的默认行为请加 `--all`。

### Step 1: 提取元数据 (extract)

从 GenBank 文件中提取记录与标本的元数据（不提取序列）。标本信息只读 `source` 特征，不读 gene、rRNA、CDS 等特征。

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
| 线粒体基因 | coi, 16s, 12s, cob, mtgenome, cox3, cox2, nd1, nd2, nd3, nd4, nd4l, nd5, nd6, atp6, atp8 |
| 核基因 | 18-28s, its1-its2, ef-1, 28s, 18s, h3 |
| 叶绿体基因（植物） | plastome, rbcl, matk, trnl-trnf, psba-trnh, rpl16, rps16, ndhf, rpoc1, ycf1, rpl32-trnl, trns-trng, atpf-atph, petb-petd, psbe-petl, rbcl-atpb, trnk, rps4, atpb |

Definition 同时写有 ITS 和 18S 或 28S 的跨区段 rRNA 记录（如 "small subunit ribosomal RNA gene and internal transcribed spacer 1"、"internal transcribed spacer 2 and large subunit ribosomal RNA gene"）得到它覆盖的每个类型，如 `gene_type` 为 `18s,its1-its2` 或 `28s,its1-its2`，`Assignment_reason` 注明 "multi-region rDNA"。完整的 18S–28S 记录仍为 `18-28s`。

匹配到多个基因的线粒体记录，只有标题写明 genome（"mitochondrion, complete genome"、"mitochondrial genome" 等；单独的 "complete sequence" 不算）才归为 `mtgenome`，否则保留各基因并标记冲突。环状（circular）记录归为其细胞器的基因组（叶绿体 → `plastome`，否则 `mtgenome`）。植物（Viridiplantae 或叶绿体记录）第二轮匹配到其它细胞器的基因时丢弃（叶绿体 16S 不归入线粒体 `16s`，线粒体 rrn26 不归入 `28s`）；动物保留该基因并标记。

**输出文件**: `assigned_genes_types_all.csv`

---

### Step 2b: 按序列核对 rRNA 类型（默认）

v0.05 起，`g2t` 在第 2 步之后用 profile HMM（pyhmmer）逐条核对 rRNA、未匹配和有冲突的记录：序列结果与标题或注释不同时，改用序列的类型。

- 例：Hesionidae 的 MG523356–60 标题写 "small subunit … mitochondrial"，实际是 16S，从 12S 列移到 16S 列。
- 原来按注释的结果保存在 `labeled_genes/*.annotation.csv`；逐条结果在 `verified_rrna/verify.csv`。
- `--no-verify` 恢复只按注释归类（v0.04 的做法）。

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

不含数字的占位值（`n/d`、`s.n.`、`no voucher`、`NA`、`unvouchered` 等）视为空，改用下一个字段。

**输出文件**: `updated_species_voucher.csv`

---

### Step 3b: 凭证号核对 (reconcile)

同一标本的不同基因，提交者常把凭证号写成不同形式（`COI_ZMMU_MSU_WS399` / `28S_ZMMU_WS399` / `WS399`），Step 3 按字面匹配会拆成多行。本步骤依据记录中的证据合并：

| 证据 | 内容 |
|---|---|
| 候选 | 共享核心编号（取自每个编号字段：specimen_voucher、isolate、culture_collection、clone、strain），并且同一物种，或同一属且编号完全相同（忽略大小写和标点，如 `SIO:BIC:A6302` = `SIO-BIC A6302`；同一标本以 *Bathykurila* sp. 和 *Bathykurila* sp. AH-2023 两个名称提交）。核心编号：去掉基因名、末尾注释如 "(holotype)" 和标点；`2014-1234` 这类年份–序号整体保留；`WS399, WS400` 这类列表每个编号各算一个；≥5 位的数字去掉前缀后也可配对（`USNM 123456` ~ `123456`）。其他配对（不同物种、或物种名缺失）只报告、不合并；不同属的完全相同编号在两行的 `match_evidence` 中标注（"same voucher as … (different genus, …)"） |
| 强 | 编号完全相同（忽略大小写和标点；两组含同一基因时只算中等）、同一篇论文（不含 Direct Submission）、完整采集日期相同、坐标相差 ≤0.01° |
| 中 | 第一作者相同、采集人相同、详细地点相同、坐标相差 ≤0.5°、两组含同一基因时编号完全相同 |
| 矛盾（阻止合并） | 两组都有 specimen_voucher 且毫无共同编号（即使 isolate 相同，如采集事件号）；共享核心后接不同的尾号（批号或子样本号：`BIOUG06836-E01` / `-E02`、`SIO:BIC:A9392-1` / `-2`），或只有一边多出数字尾号（`CE438` / `CE438_2`；多出的是含字母的编号如 `WS16799` / `WS16799_Pr21` 不算矛盾）；共享核心两边是不同的机构、采集人或样品编号（`USNM 123456` / `MCZ 123456`、`Juslen s.n. 2003` / `Muola s.n. 2003`、`Polychaeta1D9_St95` / `Polychaeta1G8_St95`；一边没有这类词时照旧可合并：`WS399` / `COI_ZMMU_MSU_WS399`、`H & W 16149` / `C.B. Hellquist, J.H. Wiersema 16149`）、日期不符、国家不同、坐标相差 >0.5° |

判定：有强证据 → 合并（high）；中等证据 ≥2 → 合并（medium）；否则不合并。两组已含同一基因时，只有采集日期相同或坐标相同才能合并（同一篇论文不够，因为一篇论文报道许多标本）；这样未合并的配对在 `match_evidence` 中标注（"similar voucher … check"）。合并两个标本前会检查它们的每一对记录，所以 A–B–C 这样的链条不会把互相矛盾、物种不同（上述同属规则除外）、或会让同一基因在一个标本中出现两次而又没有日期/坐标证据的 A 和 C 并在一起。

一个凭证号可能是一个批次或一张标本台纸：同一凭证号下不同 isolate 测了同一基因（MCZ 70173 的 7 条 COI）时，每个 isolate 一行，`match_evidence` 注明 "one voucher, several individuals"。合并后的新名字取这些记录中出现最多的物种名。

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

| species_voucher_new | organism | organism_names | coi | 16s | 18s | ... |
|---------------------|----------|----------------|-----|-----|-----|-----|
| SpeciesA_voucher1 | Species A | | LC123456 | LC123457 | LC123458 | |
| SpeciesB_voucher2 | Species B | Species B; Species B sp. 1 | LC123459 | | LC123460 | |

`organism` 取该标本记录中出现最多的名字；`organism` 之后的 `organism_names` 列在名字不一致时列出全部名字，否则为空。

单元格内为该标本对应基因的 LocusID，多个序列用分号分隔。

**输出文件**: `organized_species_voucher.csv`

**按科（或目、属）拆分**：`g2t ... --split_by Family --drop_empty_gene_columns` 另在 `organized_genes/by_Family/` 下每科写一个矩阵（列与总矩阵相同），
`_summary.csv` 列出每科的标本数、物种数、≥2 个基因的标本数和各基因标本数；没有科名的行放在 `unassigned_Family.csv`。

---

### Step 5: 污染和鉴定错误提示（默认）

`g2t` 生成矩阵后，对每个标记（COI、16S、18S……，也包括线粒体基因组中的基因）做快速的距离检查（k-mer，不需要外部软件），提示以下几类序列：

| 提示 | 含义 | 例子 |
|---|---|---|
| `near_other_genus` | 与另一个属几乎相同（≤ 2%），而离本属较远。只用于 COI、16S 等快速进化的标记 | OQ417162 "*Polyophthalmus pictus*" 与 *Armandia* 相差 0.2% |
| `same_species_far` | 同名的两条记录相差超过 15%。提示中写明各自离本属多远，离本属远的那条更可疑；也可能是隐存种 | HM746718 "*Ophelia neglecta*" 16S 与同种线粒体基因组 OZ260669 相差 24.3% |
| `ncbi_unverified` | NCBI 自己标注无法核实：标题以 "UNVERIFIED:" 开头，或 COMMENT 写 "unable to verify source organism / sequence" | 海女虫科 28 条 |
| `remote` | 与同一标记的所有序列都相距很远（不相关的基因，或混入了远缘类群） | |

结果写入 `organized_genes/screen.csv`，矩阵新增 `screen_flags` 列；`--no-screen` 关闭这一步。也可以单独运行：

```bash
g2t-screen -i gb/markers --g2t OUT
```

**建树检查（`--tree`，可选，更准）。** 距离检查只能在下载的这批序列内部比较：混入的序列若来自数据中没有的类群（例如一条名为 *Ophelia* 的序列实际来自鳞腹沙蚕科），或者序列间相差已超过 20%，就判断不了。`--tree` 做法如下：

1. 自动从 NCBI 分类中选出近缘的几个科作为背景并下载，例如海蛹科 → Scolecida 下的 6 个科；
2. 每个标记先统一序列方向，再比对（FAMSA）、修剪、建树（FastTree），做法同 OrthoFinder；
3. 看每条序列在树上与谁聚在一起：落在别的科或别的属里、且支持率 ≥ 0.9 的序列给出提示；
4. 再用距离复核排除假象，并检查同一标本的其他基因是否位置正常。

```bash
conda install -c bioconda famsa fasttree      # 需要 FAMSA 或 MAFFT，加上 FastTree 或 IQ-TREE
g2t-screen -i gb/markers --g2t OUT --tree
g2t-screen -i gb/markers --g2t OUT --tree --background "Capitellidae,Scalibregmatidae"   # 自己指定背景科
g2t-screen -i gb/markers --g2t OUT --tree --iqtree                                       # IQ-TREE + bootstrap，更慢
```

输出在 `OUT/organized_genes/phylo/`：

| 文件 | 内容 |
|---|---|
| `phylo_flags.csv` | 记录错误：单条序列落在别的科或属中。`evidence` 列给出跨基因印证——同一标本其他基因都正常，说明这条序列本身有问题（污染或混样）；其他基因也错位，说明可能是标本鉴定错误 |
| `phylo_taxon_notes.csv` | 分类问题：整个属的位置有问题、两个属互相嵌套（异名）、本属只有这一条序列 |
| `<基因>.tre` | 每个标记的树 |

背景序列下载到 `OUT/background/`，再次运行时只下载新增的记录。

海蛹科实测（约 3 分钟，含下载）：只提示 PX620130（落在 *Hyboscolex* 里）、HM746718（落在 *Notomastus* 里）、OQ417162（落在 *Armandia* 里），三条都与人工核查一致。

提示只表示值得核查，请对照原文、标本和系统树判断。评估见 `docs/evaluation/screen.md` 和 `docs/evaluation/phylo-screen.md`。

### 矩阵中的凭证号

矩阵在 `organism`、`organism_names` 之后有三列：`voucher_as_submitted`（该标本在 GenBank 中出现过的全部 specimen_voucher 原样写法，用于看出提交差异）、`voucher_standardized`（每个标本一个统一凭证号，只由 specimen_voucher 生成，只有 isolate 等编号的标本此列为空；去基因名前缀、取信息最全的写法，按 INSDC `/specimen_voucher` 格式 `[机构代码:[收藏代码:]]标本号` 书写，如 `COI_ZMMU_MSU_WS2585; 28S_ZMMU_WS2585; WS2585` → `ZMMU:MSU:WS2585`，`ZMMU:WS30980; ZMMU_WS30980` → `ZMMU:WS30980`）、`voucher_note`（几种写法不只是前缀/分隔符不同时提示人工核对）。

### Step 6: 校正元数据 (curate，可选)

#### 目的

序列通常在投稿时提交到 GenBank。审稿后常有改动：标本重新鉴定或组合变更、被描述为新种、题目改变、坐标和日期更正。GenBank 记录只能由提交者更新。INSDC 欢迎"作者对错误的更正和记录的更新"（[INSDC policy](https://www.insdc.org/policy/)）；NCBI 要求文章发表后由提交者把出版信息发给 GenBank（[GenBank 概述](https://www.ncbi.nlm.nih.gov/genbank/)），来源信息按表格提交更新（[Update GenBank records](https://www.ncbi.nlm.nih.gov/genbank/update/)）。实际上很多记录从未更新，下载的数据集里元数据可能过时或错误。公共序列鉴定错误和注释不全已有文献记载（Bridge et al. 2003, *New Phytologist* 160: 43–48, [doi:10.1046/j.1469-8137.2003.00861.x](https://doi.org/10.1046/j.1469-8137.2003.00861.x)；Meiklejohn et al. 2019, *PLOS ONE* 14: e0217084, [doi:10.1371/journal.pone.0217084](https://doi.org/10.1371/journal.pone.0217084)）。

Step 6 让分类学者把核实过的信息用到"标本 × 基因"矩阵上，不修改 GenBank，也不丢失原始值：

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

## 线粒体基因组统计（--header-only、g2t-mitostats、g2t-cite）

只想知道一个类群的线粒体基因组用什么测序、什么软件组装、发表在哪里，不需要下载序列：

```bash
g2t-download -t Annelida -o gb --mitogenome --header-only   # 只下记录头，约为完整记录的 1/12
g2t-cite -i gb/mitogenome-header -o mito/cite.xlsx          # 可选：到 Crossref 找 "Unpublished" 记录的论文
g2t-mitostats -i gb/mitogenome-header -o mito --citations mito/cite.xlsx
# 打开 mito/mitostats_report.html
```

- 输出：`mitostats.csv`（每条记录一行，原文与归类并排）、按年 / 组装软件 / 期刊的计数表、累计面积图和冲积图（PNG/PDF）、单文件离线交互报告。
- 规则（`src/g2t/data/mito_methods.yaml`）：有三代平台记为 Long-read（含混合组装），否则有二代平台记为 NGS，否则 Sanger；**没有 Assembly-Data 注释的记为 Unknown，不当作 Sanger**。年份 = 发表年，否则提交年，否则建档年（`year_source` 列）。
- `g2t-cite` 分级：high / medium 进入 `cite.xlsx` 的 *Curate* 页（可直接给 `g2t-curate -u` 和 `g2t-mitostats --citations`），check 行在 *Candidates* 页标黄，需人工确认。
- 重跑只下载新记录、只查询新的论文组。
- 完整实例（环节动物 745 条，含每一步的真实输出）：[docs/tutorials/mitogenome-stats.md](docs/tutorials/mitogenome-stats.md)；全后生动物及增量更新：[docs/tutorials/metazoa.md](docs/tutorials/metazoa.md)。

## 用序列核验基因归属（g2t-verify）

g2t 按记录标题判断基因类型，标题写错就会归错。例如 Hesionidae 有 5 条 16S 标题写作 "small subunit ribosomal RNA; mitochondrial"，被归到 12S。`g2t-verify` 用 12S、16S、18S、5.8S、28S 的 profile HMM 扫描序列本身，定位各区段，并由 5.8S 推出 ITS1、ITS2：

```bash
pip install pyhmmer
g2t-verify -i gb/markers --g2t out -o verify
```

- **结论**：agree 一致；partly 部分一致（例如标题含 28S，但序列只有 ITS）；different 是别的 rRNA；no_hit 标题说是 rRNA，但没有模型命中；unverifiable 只有间隔区，无法定位；new 是 g2t 未归类的 rRNA。
- **输出**：`verify.csv`；`regions/*.fa`，即按坐标切好的各区段序列。
- **数据和评估**：见 [docs/evaluation/rdna-regions.md](docs/evaluation/rdna-regions.md)。

## 输出文件说明

```
output/
├── gb_metadata/                          # Step 1 提取
│   ├── final.csv                         # 主输出：元数据 + 组装信息，每条记录一行
│   ├── metadata.csv                      # source 限定词、文献、分类信息（final.csv 的一部分）
│   ├── assembly.csv                      # 组装 / 测序注释（final.csv 的一部分）
│   ├── extraction_report.json            # 读入、失败、无效的文件与记录数；错误信息
│   ├── record_tracking.csv               # 仅 g2t-extract：每条解析的记录一行（文件、序号、登录号、状态）
│   ├── record_errors.csv                 # 仅 g2t-extract，有时才有：无法解析的记录
│   ├── taxonomy_cache.json               # TaxID → 纲 / 目 / 科 / 属（之后的运行复用）
│   └── failed_extract/                   # 无法读取的输入文件副本（有时才有）
├── labeled_genes/                        # Step 2 基因类型
│   ├── assigned_genes_types_all.csv      # 主输出：有基因类型的记录（第 1 + 2 轮）
│   ├── assigned_genes_type.csv           # 第 1 轮（标题关键词、长度范围）
│   ├── assigned_genes_types2.csv         # 第 2 轮（未匹配记录按所有区室重新检查）
│   ├── multiple_matched.csv              # 一条记录对应多个基因类型（如 18s,its1-its2）
│   ├── conflicted_genes.csv              # 标题与基因类型矛盾（如 "18S ...; mitochondrial"）
│   ├── unmatched_sequences.csv           # 第 1 轮后无基因类型
│   ├── unmatched_sequences2.csv          # 第 2 轮后无基因类型
│   ├── unmatched_report.txt              # 未匹配记录按类别统计，附示例
│   ├── unmatched_categorized.csv         # 未匹配记录及其类别
│   ├── filtered_records.csv              # 被过滤掉的记录及原因（长度、缺长度、重复 LocusID 等）
│   └── record_status.csv                 # 每条输入记录一行：assigned / unmatched / filtered + 原因
├── updated_species_vouchers/             # Step 3 标本；Step 3b 核对
│   ├── updated_species_voucher.csv       # Step 3 输出：带标本编号的记录
│   ├── reconciled_species_voucher.csv    # Step 3b 输出：凭证号核对后（match_* 列）
│   ├── reconcile_report.csv              # Step 3b：每对候选的证据与判定
│   ├── shorten_species_voucher.csv       # 只含编号相关列，便于检查
│   ├── only_haplotype.csv                # 带单倍型名称的记录
│   ├── species_voucher_groups.csv        # 被 2–8 条记录共用的编号，含 frequency_group 列
│   ├── species_voucher_new_group_N.csv   # 同上，每个组大小 N 一个文件
│   └── count_species_voucher.txt         # 编号出现 1、2、3 …… 次的个数
├── organized_genes/                      # Step 4 矩阵
│   └── organized_species_voucher.csv     # 主输出：每个标本一行，每个基因一列
├── .g2t_params.json                      # 上次运行的核对设置（设置改变时 --resume 重跑 3b / 4）
└── g2t_run.json                          # 完整流程：日期、命令、各软件版本、全部参数、输入 GenBank 文件（路径、大小、SHA-256）、分类信息来源
```

---

## 常见问题

### Q: 如何处理没有 specimen_voucher 的序列？

程序会依次尝试 isolate → culture_collection → clone → strain，最终使用 ACCESSION 作为兜底。

### Q: 如何添加自定义基因类型？

编辑 `src/g2t/data/gene_dict.yaml`（同义词的唯一来源），添加新的基因类型和同义词：

```yaml
gene_types:
  new_gene:
    synonyms:
      - synonym 1
      - synonym 2
```

或使用外部词典：

```bash
g2t-classify -i input.csv -o output --path_dict custom_dict.txt
```

`--path_dict` 接受 `.txt`（每行一个基因：`new_gene: synonym 1, synonym 2`）、`.csv`（列 `gene,synonyms`，同义词用逗号分隔）或 `.py`（定义变量 `gene_dict`）；不接受 `.yaml`。其中的条目会加入内置词典。

新基因类型还需列入 `g2t-classify --mtgenes_list` 或 `--ntgenes_list`（否则只在第二轮匹配，没有该基因自己的长度范围），并列入 `g2t-organize --gene_order`（否则矩阵中没有这一列）。`g2t` 流程命令不传递这些参数，因此需单独运行 Step 2 和 Step 4。

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
pip install pandas biopython openpyxl pyyaml
```

---

## 性能建议

- **大文件处理**: 使用 `--stream` 参数，避免内存溢出
- **多文件处理**: 使用 `--batch` 参数并行处理
- **断点续传**: 使用 `--resume` 参数跳过已完成步骤；凭证号核对的设置与上次不同时（记录在 `OUT/.g2t_params.json`），Step 3b 和 4 会重跑。同一目录重跑 `g2t-classify` 会先清除该步骤上次的输出。

---

## 不足与已知限制

1. **基因类型按 DEFINITION 整条判定。** 跨区段的 rRNA 记录（如"18S … ITS1""ITS2 … 28S"）放入它覆盖的每一列（`18s,its1-its2`、`28s,its1-its2`）：同一登录号出现在两列，代表整条记录，而不是单独的区段。匹配到其他基因组合的记录得到全部匹配类型（如 `16s,coi`、`28s,h3`），并标记 `Conflict=True`（`Conflict_reason` 为 "Multiple candidates…"）；匹配到两个以上线粒体基因且长度 ≥3,000 bp 的记录，只有标题写明 genome 时才归为 `mtgenome`。详见 [BUGS.md](BUGS.md)。
2. **不切分序列。** 线粒体基因组（`mtgenome`）的登录号会被复制到全部线粒体基因列（coi、16s、12s、cob、cox2、cox3、nd1–nd6、nd4l、atp6、atp8），叶绿体基因组（`plastome`）复制到全部叶绿体基因列，18S–ITS–28S 记录（`18-28s`）复制到 18s、28s、its1-its2，**不论该记录是否真的注释了这个基因**（如 FN689349 没有 12S 注释，仍出现在 12s 列）。可用 `g2t-organize --mtgenome_includes coi 16s` 限定列。不会切出各基因的序列；比对前需另行提取。
3. **长度过滤。** 默认长度 150–50,000 bp 以外（以及缺长度、LocusID 重复）的记录在分类前被去掉。它们不会静默消失：逐条列在 `filtered_records.csv`（含原因），`record_status.csv` 给每条输入记录恰好一个状态。若要让它们也参与分类，三个范围都要放宽：`g2t-classify … --length_range2_all 1:1000000 --length2_mtgenes 1:1000000 --length2_ntgenes 1:1000000`。
4. **只识别 41 类基因，按关键词匹配。** 其他位点（核蛋白编码基因、微卫星、Hox 基因等）进入 `unmatched_sequences.csv`，除非加入基因词典（见常见问题"如何添加自定义基因类型"）。DEFINITION 中的标本编号（`isolate CO2`、`voucher COI-12`）和英文单词 "its" 不参与匹配，并排除了几种形似的情况（16S rRNA methyltransferase、histone H3 lysine/demethylase、elongation factor-1 beta/gamma），但写法特殊时关键词匹配仍可能出错。
5. **凭证号与元数据。** Step 3 把同一物种中凭证号字符串完全相同的记录直接合并，不做进一步检查；不同标本恰好同号时会被误并。Step 3b 还会合并同一属中以不同名称提交、编号完全相同（忽略大小写和标点，任一编号字段）的组（*Bathykurila* sp. / *Bathykurila* sp. AH-2023）；两组含同一基因时需采集日期或坐标相同；不同属的相同编号不合并，并在 `match_evidence` 中标注。Step 3b 依赖提交到 GenBank 的元数据：没有共同的论文、日期、坐标或采集人时，同一凭证号的不同写法不会合并（宁缺毋滥）；两组已含同一基因时，同一篇论文不够，该配对改为标注。只差批号或子样本尾号的编号（`CE438` / `CE438_2`），或编号两边的机构、采集人、样品编号不同（`USNM 123456` / `MCZ 123456`），视为不同标本。矩阵中的地点、日期等元数据逐列取该标本第一条有值的记录，因此可能来自不同记录；任一记录被标记冲突，`Conflict` 即为 True。
6. **物种名按提交原样。** 不与 WoRMS、NCBI 异名等分类权威核对；错误鉴定或过时名称保持原样，除非在 Step 6 中校正。
7. **下载。** 只查 NCBI nuccore（不含 BOLD、仅在 ENA 的数据、SRA）。`--gene` 按基因字段和标题关键词检索，写法特殊的记录可能漏检；按基因检索也会返回含该基因的线粒体基因组。
8. **格式损坏的文件。** GenBank 文件中某条记录无法被 Biopython 解析时，该文件中其后的记录不会被读取；`extraction_report.json` 中该文件标为部分完成，并在警告中给出应有与实际解析的记录数。
9. **校正（Step 6）。** B 只校正 A 中已有的标本。没有 GenBank 记录的标本列为"匹配不到"，不会加入 C。空单元格表示"不改"，因此不能通过 B 清空 A 中的值。列识别是根据列名和"像登录号的值"的启发式判断，请查看 *Column mapping* 表，识别错时用 `--map` 指定。B 中的值被视为正确：程序不会把 B 本身与文献或 WoRMS 核对。
10. **代码状态。** v0.04 为早期版本，610 项单元测试在 Python 3.12、3.13 和 3.9 上通过（本地测试；3.9 上其中 3 项 Excel 测试需要 openpyxl）；`ruff` 与 `mypy` 无警告。

---

## 引用

g2t 目前没有发表论文，请引用软件本身及所用版本：

> Yang, D. (2026). *g2t: GenBank to Taxonomy* (version v0.04) [Computer software]. GitHub. https://github.com/deyuanyang92-dev/gb2taxonomy

```bibtex
@software{yang_g2t_2026,
  author  = {Yang, Deyuan},
  title   = {g2t: GenBank to Taxonomy},
  version = {v0.04},
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
