# 更新日志 / Changelog

## v0.02（2026-10-09）

包内版本 `0.0.2`。本版来自一次全面代码审查（新代码 19 项、原有模块 25 项，均写出复现并实际运行确认），修复后 344 项测试通过。
用 Priapulidae 数据回归：分类结果（496 条）、未匹配（22 条）、凭证号步骤输出与 v0.01 逐列一致；凭证号核对仍合并 94 个标本。

### 新增（自 v0.01 起）
- **被过滤记录不再静默丢弃**（classify）：
  - `filtered_records.csv`：每条被过滤的记录及原因——全局长度过滤（默认 150–50,000 bp）、长度缺失、`--moleculetype/--mol_type/--length/--organelle` 过滤、重复 LocusID（保留首条）。
  - `record_status.csv`：每条输入记录恰好一行，状态为 assigned（附 gene_type）/ unmatched / filtered，均附原因；第 2 轮 recheck 因长度被跳过的未匹配记录在原因中注明。
  - Priapulidae 实测：548 = assigned 496 + unmatched 22 + filtered 30（此前 30 条短于 150 bp 的 FOXA3/Hox 片段在任何输出中都看不到）。

### 修复：下载（g2t-download）
- 记录集变化后续传会丢新记录、重复旧记录：现在每次重新取 accession 列表，批文件按其中的 accession 逐条核对后才复用。
- 不同的 `--query` 或 `--include-*` 共用一个子目录并互相复用批文件：子目录名现在包含这些选项（`custom-<检索式哈希>`、`incl-wgs`、`all-…`）；同一目录换检索式时报错。
- 批大小变大或记录变少后，多余的旧批文件不再残留。
- `--dry-run` 不再创建目录或改写 manifest；manifest 在下载结束后写入，并记录 accession 列表校验值与是否完整。

### 修复：凭证号核对（g2t-reconcile）
- 同一物种同一核心编号下 ≥3 个不同标本时新键冲突、被并成一个。
- 经第三条记录的传递合并可绕过"同基因须强证据"和冲突检查：合并前检查两个标本的每一对记录。
- 跨 180° 经线的坐标被算作相距 360°。
- 核心编号：保留年份–序号，去掉末尾括号与标点，列表中的每个编号都作为候选，≥5 位数字可去前缀配对。
- 缺 organism 列时报错（以前当作同一物种）；organism 为空的记录不合并。
- 多基因类型（如 `coi,16s`）的重叠判断。
- 凭证字段优先级与 voucher 步骤一致（clone 先于 strain）。

### 修复：分类（classify）
- 只要有一行缺长度，其余行长度被放大 10 倍（"550.0" → 5500），导致误判 18-28s / mtgenome。
- 关键词误配：DEFINITION 中的标本编号（`isolate CO2`、`voucher COI-12`）和英文单词 "its" 不再参与匹配；排除 16S rRNA methyltransferase、histone H3 lysine/demethylase、elongation factor-1 beta/gamma。
- 冲突解析的结果随 Python 哈希种子变化（不可重复）。
- 同一目录重跑时读到上次的旧结果；关闭第 2 轮时不生成 `assigned_genes_types_all.csv`。
- 所有记录都被过滤时报 `KeyError` 失败；带千分位逗号的长度被当作缺失；含空格的列名被改写；以标点开头的同义词导致正则错误；`none`/`all` 长度范围导致崩溃；mtgenome 下限写死为 3000、忽略 `--length_range2_mtgenome`。
- `g2t-classify` 命令行也生成 `record_status.csv`。

### 修复：提取（extract）
- `--batch` 并行模式完全不可用（局部函数无法序列化），现可用并输出 `final.csv`。
- 文件中途解析失败仍报成功、后续记录静默丢失：现在与 `//` 计数对账，写入 `extraction_report.json` 并警告。
- 合并 `final.csv` 时数值样式的列被改写（`007` → `7`，PubMed `32671913` → `32671913.0`）。
- Python 3.13 下 ete3 无法导入，Class/Order/Family/Genus 全空：改为从 NCBI Taxonomy 获取并缓存。
- 非 taxon 的 db_xref（如 BOLD）被丢弃；stream 与非 stream 输出列不一致；无 source 时 organism 为空；死计数 `assembly_failed`。ete3 不可用时不再用 `input()` 询问。

### 修复：凭证号 / 矩阵（voucher, organize）
- `--mtgenome_includes` 无效；`--normalize_columns` 把列名转小写导致第 4 步失败；`g2t-voucher` 命令行默认值与 API 不一致；纯空白/标点的凭证号不回退；`Conflict` 只取第一条记录；不在 gene_order 中的基因被静默丢弃；无标本键的记录全部并成一行 UNKNOWN。

### 修复：其他
- `g2t.reconcile()`、`g2t.organize()` 等包级函数调用一次后被同名子模块覆盖（README 的 Python API 示例因此出错）。
- `--resume` 不随凭证号核对设置的变化重跑。

### 未修复（已记入 BUGS.md）
- 第 2 轮复查中 Topology=circular 直接判为 mtgenome，且忽略 `--mtgenes_list/--ntgenes_list`（L-18）。
- `g2t` 命令的 `--extract_extra/--classify_extra/--voucher_extra/--organize_extra` 解析后未传给流水线，目前无效。
- 格式损坏的记录之后同一文件的记录不读取（已报告，未恢复）。

### 文档
- README/USAGE：下载续传与目录命名、旧下载脚本的差异、放宽长度需同时设三个范围、mtgenome 的列复制规则与 `--mtgenome_includes`、矩阵元数据逐列取值、格式损坏文件的处理。

## v0.01（2026-10-08）

首个版本号。包内版本 `0.0.1`（PEP 440 会把 "0.01" 规范化为 "0.1"，故包内写 0.0.1；GitHub 标签为 `v0.01`）。

### 新增

**Step 0 下载 `g2t-download`（`g2t/download.py`，`g2t/ncbi_genes.py`）**
- 按类群（名称或 taxid）从 NCBI nuccore 下载 GenBank 记录。
- 默认 "markers"：排除 WGS contig、mRNA、RefSeq（预测模型及 NC_ 等与 INSDC 重复的拷贝）。实测 Priapulidae 共 95,726 条，其中 WGS 69,704、mRNA 23,053；标记基因集 548 条。`--all` 可下全部。
- 选择：`--mito`、`--mitogenome`（10–30 kb）、`--gene COI,18S,…`（18 个基因，见 ncbi_genes.py）、`--minlen/--maxlen`、`--query`、`--include-wgs/--include-mrna/--include-refseq`。
- `--dry-run` 只报告类群构成与各基因条数。
- 按 accession 列表分批下载（不依赖会过期的 WebEnv），每批校验条数，指数退避重试，重跑即续传。
- 每种选择独立子目录 `<out>/<tag>/`，附 `accessions.tsv` 与 `manifest.json`（检索式、日期、条数）。
- `scripts/download_entrez.py` 改为该模块的包装，旧参数兼容。

**Step 3b 凭证号核对 `g2t-reconcile`（`g2t/reconcile.py`）**
- 问题：同一标本不同基因的凭证号常写法不同（`COI_ZMMU_MSU_WS399` / `28S_ZMMU_WS399` / `WS399`），Step 3 按字面匹配拆成多行。
- 候选：同一物种 + 核心编号相同（去掉基因名与馆藏缩写，≥4 字符）；不同物种同号只报告、不合并。
- 证据：强 = 同一论文（不含 Direct Submission）/ 同采集日期 / 坐标 ≤0.01°；中 = 同第一作者 / 采集人 / 详细地点 / 坐标 ≤0.5°；矛盾 = 日期不符 / 国家不同 / 坐标 >0.5°。
- 判定：有矛盾不合并；≥1 强 → 合并 high；≥2 中 → 合并 medium；两组已含同一基因时只接受强证据。
- 输出 `reconciled_species_voucher.csv`（新增 species_voucher_g2t / voucher_core / match_basis / match_confidence / match_evidence）与 `reconcile_report.csv`（每对候选的证据与判定）。
- 完整流程默认开启；`g2t --skip_reconcile` 关闭，`--reconcile_min_confidence high` 只按强证据合并。
- Priapulidae 实测：130 对候选均有"同一论文"强证据，94 个标本合并（high）；WS3020（*Halicryptus spinulosus* 28S vs *Priapulus caudatus* COI/16S）跨物种，未合并。

### 修改
- `organize`：最终矩阵增加 match_basis / match_confidence / match_evidence 三列。
- `_pipeline` / `_cli`：流程为 extract → classify → voucher → reconcile → organize；新增 `g2t-download`、`g2t-reconcile` 命令。

### 测试
- 新增 38 项（test_reconcile.py 22、test_download.py 16），共 232 项通过；新增文件通过 ruff 与 mypy。

### 未改动（有意保留）
- **18S / 28S / ITS 仍按记录的 DEFINITION 整条归类**（原逻辑）。已知限制：跨区段记录归入 its1-its2，例如 AY210840、AH010828（5.8S+ITS2+约 3.66 kb 28S）以及 18S+ITS1 的记录（1.4–1.8 kb，主体为 18S）。拆分方案（先按特征表坐标切，切不开再用 ITSx）已评估，暂不实施。见 BUGS.md。
