# 更新日志 / Changelog

## 未发布 / Unreleased

### 下载（g2t-download）
- **默认跳过 >100 kb 的核基因组记录**（染色体级组装、基因组 scaffold、TSA 主记录、大片段克隆；`--include-large` 保留，显式 `--minlen/--maxlen` 时不加此限制）。**细胞器记录不受长度限制**（`mitochondrion[filter]`、`chloroplast[filter]` 或标题含 mitochondrion/mitochondrial/chloroplast/plastid）：环节动物线粒体基因组最大约 47.5 kb（*Glycera fallax* OZ210882），植物常 >100 kb。Polynoidae：248 条 “genome assembly, chromosome: N”（单条约 100 MB）使第 1 批达 5.8 GB、约 1 小时；跳过后每批约 6 秒，全类群约 3 分钟。预检报告新增 “nuclear > 100 kb” 条数。
- **`--mitogenome` 修正**：原检索式写死 10–30 kb 且只认 “complete genome/mitogenome” 标题，漏掉 35–48 kb 的线粒体基因组和 Darwin Tree of Life 的 “genome assembly, organelle: mitochondrion”。现在 ≥10 kb 无上限，并接受 “mitochondrial genome”/“genome assembly” 标题。Annelida：634 → 733 条（新增 99 条均为 DToL 线粒体基因组，10.3–47.5 kb）；Polynoidae：18 → 31 条。
- **失败批次自动二分**：一批反复失败（如某条超大记录导致传输中断）时拆半重试，直到定位到单条失败记录；其余记录照常入库，manifest 标记未完成。
- **按记录增量下载**：新模块 `g2t/recstore.py`，记录按 `accession.version` 存入 SQLite（默认 `<out>/_records.sqlite`，`--store` 可多个类群共用）。只下载库中没有的记录；换选择、重叠类群直接复用；新版本（`.2`）才重新下载；旧版批文件自动导入。`<tag>/batch_NNNN.gb` 由库重写，下游不变。此前 NCBI 列表头部新增几条记录，其后所有批都要重下。
- **按 NCBI 修改日期增量下载 `--since auto|YYYY-MM-DD`**：只检索上次下载（manifest 日期，往前留 3 天余量）以来新增或修改的记录（Entrez `[MDAT]`），不再取完整登录号清单；更新版本替换旧版本。局限：查不到被撤下的记录、以及记录未改动而 NCBI 分类树变化带来的增减——不加 `--since` 跑一次即完整核对（manifest 记 `last_full_check`）。Polynoidae 实测：把上次日期设为 2026-09-01，`--since auto` 只查 46 条修改记录，21 s；结果与完整清单逐条一致（5858 条、无重复）。
- **不用本地库 `--no-store`**：以选择目录中的批文件为唯一副本，扫描已有记录，只下载缺失的写成新批文件，旧版本/撤下的记录从旧文件中删除。可与 `--since` 组合。Polynoidae 从零 56 s。
- **变化记录**：与上次列表比较，新增/更新/撤下写入 `changes.tsv`，并追加到 `changes_history.tsv`；g2t 默认检索式变化时标注 "(query changed)"。
- **传输提速**：efetch 改为 gzip 传输（实测 500 条 2.26 MB → 0.29 MB），60 s 无数据即放弃重试（原 Bio.Entrez 无超时，断线要等很久）；每批默认 500 条（实测与 200 条耗时相近，吞吐约 2.5 倍）。Polynoidae 全部 5858 条从零下载 138 s、20 MB。
- **正式下载不再做统计检索**：类群构成与各基因条数需约 23 次 esearch（大类群每次 2–15 s；Nereididae 合计 86–380 s），下载本身用不到，现只在 `--dry-run` 或 `--report` 时做；下载后各基因条数由 pipeline 统计。所有 Entrez 请求加 60 s 超时。Nereididae 17,403 条全复用重跑 127 s → 37–100 s（波动来自 NCBI/网络）。
- 并行请求 `-w/--workers`（默认 3），登录号清单分页也并行获取；`--dry-run` 报告已存/待下载条数；选择未变且无新记录时不重写批文件。
- **大规模实测**（Phyllodocida，188,739 条）：复用此前 Nereididae 的 17,403 条，新下载 171,336 条用时 835 s（约 205 条/s），途中 40 次网络错误全部自动重试，0 失败；全程 18 min，内存 744 MB，记录库 314 MB。全部复用的重跑 86 s。
- `DownloadResult` 字段改为按记录计数：`fetched`、`reused`、`to_fetch`、`new`、`updated`、`removed`、`query_changed`（原 `downloaded`/`skipped` 按批计数，已移除）。
- 只有用 `--tag` 命名的目录在检索式变化时报错；自动命名的目录按新检索式重写。
- 测试：新增 22 项，共 410 项通过。

## v0.03（2026-10-09）

包内版本 `0.0.3`。本版新增统一凭证号和元数据校正（Step 5：NCBI 矩阵 A + 你自己的表 B → 校正后的矩阵 C）。

### 新增
- **统一凭证号**（organize）：矩阵在 `organism` 之后新增 `voucher_standardized`（去基因前缀、取最完整写法，按 INSDC `/specimen_voucher` 格式 `机构代码:收藏代码:标本号` 书写，如 `ZMMU:MSU:WS2585`）、`voucher_as_submitted`（GenBank 原样，全部写法）、`voucher_note`（写法差异超出前缀/分隔符时提示）。Priapulidae：213 个有凭证号的标本，同一物种内标准化凭证号无重复，无需提示的冲突。
- **元数据校正 `g2t-curate`**（Step 5，`g2t/curate.py`）：按 accession 和/或凭证号（任何写法，可加 `organism_match`）定位标本，用校正表更新经纬度、物种名、地点、出版物等；输出校正后矩阵、修改记录、问题清单；`--template` 生成已填当前值的模板，并另存 baseline 工作表：只有改过的单元格才算修改，未改单元格不会覆盖 GenBank 之后更新的值，改过且 GenBank 也变了的列为冲突。
- **用你自己的表校正（A + B → C）**：`g2t-curate -m A.xlsx -u B.xlsx -o C.xlsx`。B 为任意格式的表（中英文列名），自动识别凭证号/标本号、多个登录号列、拉丁名、纬度+经度（合并为 GenBank `lat_lon`）、采集日期（转为 `DD-Mon-YYYY`）、题目/作者/期刊等；`--map` 强制指定；其余列以 `user:<列名>` 加入。A 可直接用 Excel（自动读 `Matrix` 工作表）。C 含 Matrix（修改处标黄）/ Changes / Problems / Column mapping / Matrix (GenBank)。

### 修改
- **凭证号核对**：复合凭证号中每个"字母+数字"编号都作为候选（`ZMMU MSU WS14906` 与 `COI_ZMMU_MSU_WS14906_XZ5507` 现在会被比较）。Priapulidae：合并标本 94 → 100（新增的 8 对均有同一论文强证据），矩阵 323 → 315 个标本。
- **安全**：写 Excel（校正结果 C、模板）时所有数据单元格按文本写入；此前以 `=` 开头的 GenBank/用户值会被写成公式，打开时执行。
- **文档**：README / USAGE 说明 Step 5 的目的与原则（GenBank 记录只能由提交者更新，发表后常未更新；引用 INSDC 政策、NCBI 更新说明、Bridge et al. 2003、Meiklejohn et al. 2019），并在"不足"中列出校正的限制。
- 测试：新增 45 项，共 389 项通过。

### 计划（后期）
- 字段采用 Darwin Core 命名；导出 NCBI 记录更新表（source modifiers）。

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
