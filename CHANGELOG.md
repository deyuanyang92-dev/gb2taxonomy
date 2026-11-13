# 更新日志 / Changelog

## 未发布

- **`g2t-screen --tree`：用基因树检查污染和鉴定错误。** 做法参考 OrthoFinder 3：先宽收，再由系统发育决定归属；FAMSA 比对、修剪、FastTree 建树，可选 IQ-TREE。
  - **背景类群：** 自动从 NCBI 分类中选取近缘科（例如海蛹科 → Scolecida 的 6 个科），并增量下载。
  - **判断方式：** 用无根树中高支持率的分割判断位置，不受外群选择影响。
  - **区分两类问题：** 记录错误（单条序列落在别处），与分类问题（整个属的位置、异名、本属只有一条序列）。
  - **复核：** 用距离复核，排除长枝吸引等假象；同一标本的其他基因也一并检查。
  - **海蛹科验证：** 只提示 PX620130（落在 *Hyboscolex* 里，相差 0.2%）、HM746718（落在 *Notomastus* 里）、OQ417162（落在 *Armandia* 里），正是人工核查确认的三条，没有误报。
  - 详见 `docs/evaluation/phylo-screen.md`。

## v0.06（2026-11-13）

包内版本 `0.0.6`。

- **第 1 步提取更多记录头信息**（新列）。这些信息原来都没有提取，而测序方法和 NCBI 的核验意见正写在这里：
  - `Comment`：自由文本 COMMENT；
  - `BioProject`、`BioSample`、`SRA`：取自 DBLINK；
  - `StructuredComment`：Assembly-Data 以外的结构化注释；
  - `primer_fwd`、`primer_rev`：从 `/PCR_primers` 解析出的引物名，多套引物用 "; " 分隔；原始的 `PCR_primers` 照旧保留。
- **`g2t-screen` 新增 `ncbi_unverified` 提示**：标题以 "UNVERIFIED:" 开头，或 COMMENT 写 "GenBank staff is unable to verify source organism / sequence and/or annotation"，即 NCBI 自己存疑的记录。
  - Hesionidae 有 28 条。
  - 按距离提示的记录中，有 6 条同时被 NCBI 标为存疑（PX640331、PX640352、PX640417、PX640477、MW278311、PP652214），这是对距离提示的独立印证。

- **`g2t-mitostats` 从更多记录头信息判断测序方法**。没有 Assembly-Data 注释时，依次看：
  1. 自由文本 COMMENT 中写明的平台。例如 Darwin Tree of Life、ATLASea 等基因组计划写的 "based on 91x PacBio data and Arima2 Hi-C data"，判为 Long-read。"Wellcome Sanger Institute" 是机构名，不当作 Sanger 测序。
  2. DBLINK 链接了 SRA 原始数据，或 COMMENT 中有 SRA 运行号（第三方组装），判为 NGS。
  3. 论文标题写明高通量方法，判为 NGS。
  4. 2012 年以前提交，推断为 Sanger。

  结果：
  - 后生动物（不含人类）的 Unknown 从 21,138 条降到 15,866 条，剩下的约 15,000 条记录头里没有任何 COMMENT；
  - Annelida 的 Long-read 从 5 条增至 95 条，Unknown 从 160 条降到 61 条。

- **`g2t-screen` 新增 `same_species_far` 提示**：同名的两条记录在快速进化的标记上相差超过 15%，就提示这两条。
  - 提示中写明各自到本属最近序列的距离，离本属远的那条更可疑。
  - 例：HM746718 "*Ophelia neglecta*" 的 16S 与同种线粒体基因组 OZ260669 相差 24.3%，最接近的是 *Ophelina*。
- **线粒体基因组参与检查**：蛋白基因取自注释；12S/16S 没有注释时，用 g2t-verify 的 HMM 定位。
- **距离计算改进**：改为按较短序列计算共享 k-mer，片段与完整基因比较时，不会再因长度不同而算大。

## v0.05（2026-11-13）

包内版本 `0.0.5`。本版重点：

1. **按序列归类 rRNA，默认开启。** 流程新增第 2b 步：12S、16S、18S、5.8S、28S、ITS 按序列本身归类，不再只看标题和注释。
2. **自动提示疑似污染或鉴定错误的序列，默认开启。** 流程新增第 5 步 `g2t-screen`。

测评见 `docs/evaluation/v0.05.md`。

- **流程第 2b 步：rRNA 按序列归类（默认）**
  - 第 2 步按注释归类后，`g2t-verify` 用 profile HMM 逐条核对 rRNA、未匹配和有冲突的记录。序列结果与注释不同的，改用序列的类型。
  - 对照：原注释保存在 `labeled_genes/*.annotation.csv`；表中新增 `annotation_type`、`sequence_type`、`verify_verdict` 三列；`Assignment_reason` 写明改动。
  - 逐条结果和切好的区段在 `verified_rrna/`。
  - `--no-verify` 恢复 v0.04 的做法，只按注释归类。
  - pyhmmer 改为必需依赖（Windows 除外，HMMER 没有 Windows 版本）。缺 pyhmmer 或模型下载失败时，跳过本步并给出提示。
  - Hesionidae：MG523356–60 从 12S 列移到 16S 列；只含 44–48 bp 18S 末端的 ITS 记录移出 18S 列。
- **流程第 5 步 / `g2t-screen`（新命令，默认运行）：提示疑似污染或鉴定错误的序列**
  - 对每个标记分别计算 k-mer 距离（正反链，不需比对，不需外部程序）。
  - 提示两类序列：与另一个属几乎相同（≤ 2%）而离本属较远的序列（`near_other_genus`，只用于 COI、16S 等快速进化的标记）；与同一标记所有序列都相距很远的序列（`remote`）。
  - 结果写入 `organized_genes/screen.csv`，矩阵新增 `screen_flags` 列。`--no-screen` 关闭。
  - Opheliidae COI 只提示 OQ417162（"*Polyophthalmus pictus*"，与 *Armandia* 相差 0.4%），即此前人工建树发现的那条。
  - Enchytraeidae 9,188 条提示 10 条；Hesionidae 0 条。
  - 见 `docs/evaluation/screen.md`。
- **下载加速**
  - 自动读取 API key：`$NCBI_API_KEY`，或 `~/.config/g2t/credentials` 中的 `api_key=`。
  - 有 key 时默认 10 个并发。环节动物 6,000 条实测：3 个并发 152–173 条/秒，8 个 254–258 条/秒，10 个 331 条/秒。
- **`g2t-mitostats` 不再把早年线粒体基因组记为 Unknown**
  - NCBI 从 2012 年才有 Assembly-Data 注释（Sequencing Technology）：150,025 条后生动物线粒体基因组中，2012 年以前提交的没有一条带这个注释。所以“没有注释”对早年记录说明不了测序方法。
  - 新规则（新增 `seq_basis` 列说明依据）：
    1. 有注释：按注释。
    2. 没有注释、论文标题写明高通量方法（454、Illumina、NGS、genome skimming、shotgun 等）：NGS。标题中最早的这类记录出现在 2007 年。
    3. 没有注释、2012 年以前首次提交：Sanger。在最早有注释的 2012 年，Sanger 195 条、NGS 66 条。
    4. 其余（2012 年及以后，没有注释）：仍为 Unknown，因为当时两种方法都常见。
  - 首次提交年份取 esummary 的建档日期，没有时取 Direct Submission 的 "Submitted (日期)"。
  - 已知局限：2006–2011 年用了 NGS、但注释和标题都没写明的记录，会被算作 Sanger。
- **修复**：运行时间偶尔显示为负数（如 "Time: -4.1s"）。WSL 下系统时钟会跳动，计时改用 `time.monotonic()`。


- **线粒体基因组统计 `g2t-mitostats`（新命令）**：汇总一个类群线粒体基因组的测序技术（Sanger / NGS / Long-read / Unknown）、组装软件、发表期刊和年份。输出逐条记录表 `mitostats.csv`（原文与归类并排）、按年 / 软件 / 期刊的计数表、论文式累计面积图和年份 → 测序 → 组装冲积图（PNG/PDF），以及单文件离线交互报告 `mitostats_report.html`（年份、测序类别、分类阶元、组装软件、发表状态筛选，图表联动，表格可排序并链接 NCBI，可下载筛选后的 CSV，深色模式，手机可用）。归类规则在 `g2t/data/mito_methods.yaml`，收录 GenBank 中实际出现的拼写（含 getorangelle、GetOrgnelle 等错误写法）；没有 Assembly-Data 注释的记录记为 Unknown，不当作 Sanger。依赖 `pip install -e ".[report]"`（matplotlib、plotly）。
- **`g2t-download --header-only`（新选项）**：只下载记录头（efetch `seq_start=1&seq_stop=1`），不含序列；真实长度、日期和 taxid 由 esummary 写入 `summary.tsv`。50 条环节动物线粒体基因组实测：0.10 MB vs 完整记录 1.24 MB。记录头存于单独的 `_records_header.sqlite`，之后的完整下载不会误用。重跑只下载、只汇总新的 accession.version；esummary 并行（`-w`），失败的部分留待下次补上。
- **`g2t-cite`（新命令）**：为 GenBank 中仍为 "Unpublished" 的记录到 Crossref 按标题和第一作者查找已发表论文，按标题相似度、作者重合度和物种名是否一致分为 high / medium / check / low。输出 `cite.xlsx`：*Curate* 页（high + medium，一行一个登录号）可直接作为 `g2t-curate -u` 的表 B 和 `g2t-mitostats --citations` 的输入；*Candidates* 页列出全部候选，check 行标黄。Crossref 原始结果缓存于 `crossref_cache.json`，重跑只查新出现的论文组。默认不发送邮箱（`--mailto` / `CROSSREF_MAILTO` 可选）。环节动物实测：224 条未发表记录、81 组 → high 16 组（124 条）、medium 9 组（11 条），人工逐条核对标题、作者和物种名均正确。
- **`g2t-mitostats --citations cite.xlsx`**：用 g2t-cite 结果补全期刊和发表年（只改 GenBank 仍为 Unpublished 的记录，`journal_source` 列标明来源）；Crossref 期刊全名自动对应到 GenBank 已有缩写（如 Mitochondrial DNA Part B → Mitochondrial DNA B Resour），单词期刊名须完全相同（Gene ≠ Genes）。环节动物发表比例 39% → 57%。
- **`g2t-download` 借鉴 kingfisher 的做法**（对照其文档：`--run-identifiers-list`、`annotate`、进度条、`NCBI_API_KEY`）：
  - `-a/--accessions FILE`：按登录号清单下载（一行一个，或 .csv/.tsv 第一列），代替 `-t`；没写版本号的取当前版本；记录库、`--header-only`、重跑增量照常。
  - `--summary-only`：不下载记录，只写 esummary 表 `summary.tsv`（登录号、长度、建档/更新日期、taxid、标题），最快的计数方式。
  - 进度行显示已下载记录数、每分钟记录数和预计剩余时间。
  - 设置了 API key 时默认并行数为 8（NCBI 允许 10 次/秒），否则仍为 3；`-w` 可覆盖。
  - 未采纳：多数据源自动切换。GenBank 记录只能经 NCBI E-utilities 获取；ENA 检索接口实测对线粒体基因组全文检索 2 分 12 秒后返回 HTTP 500，不能作为替代来源。
- **`g2t-verify`（新命令）**：用序列本身核验 rRNA 的基因归属。以 profile HMM（HMMER，经 pyhmmer 调用，`pip install pyhmmer` 即可，不需另装程序）扫描线粒体 12S、16S 和核 18S、5.8S、28S，定位各区段，由 5.8S 推出 ITS1、ITS2，逐条给出 agree / partly / different / no_hit / unverifiable / new，并按坐标切出各区段 FASTA。核模型取自 Rfam 15.1 协方差模型内嵌的 HMM（CC0，随 g2t 打包）；线粒体模型为 barrnap 0.9 `mito.hmm`（GPL-3.0，首次使用时下载并校验 SHA-256，不打包）；`--models` 可用自建模型。重叠的不同模型命中互相竞争；无命中的记录以 E ≤ 1e-3 重扫。实测：Enchytraeidae 600 条已知类型记录中 591 条一致，12S/16S/18S 无一混淆；Hesionidae 发现 5 条标题和注释都写 "small subunit … mitochondrial" 的 16S（MG523356–60，与 16S 记录 98–100% 相同），原冲突检查无法发现。近全长 18S/28S 每条约 0.5 秒；完整协方差模型（Infernal cmscan）太慢，未采用。
- **`g2t-verify` 加速**：Enchytraeidae 4,750 条在相同机器负载下从 11 分 36 秒（CPU 938 秒）降到 83 秒（CPU 162 秒），有缓存时重跑 9 秒；结论与每条记录跑全部完整模型的结果相比只有 4 条不同，都是改进。做法：
  - 相同序列只扫描一次，“模型 × 序列块”任务用满全部核心；
  - 长模型切成 300 个节点、互相重叠的分段模型（pyhmmer Builder 重建，E 值经过校准），分段没对上的末端只在那一小段上用完整模型补齐；
  - 先跑 g2t 所标类型的模型，没解释完的记录再跑其余模型；
  - 命中按序列和模型缓存（`~/.cache/g2t/verify_hits.sqlite`，`--no-cache` 关闭）；
  - E 值固定按每 1 Mb 序列计算，结果不随批次大小变化；
  - GenBank 改为整文件读入后切分，比 Biopython 快约 20 倍，输出相同。
  - 详见 `docs/evaluation/rdna-regions.md` 第 8 节。
  - 新选项：`--no-cache`、`--all-models`。
- **全后生动物实测后的修正**（150,025 条线粒体基因组记录，2026-11-13）：
  - 交互报告在 15 万条记录时打不开：`Math.min(...years)` 把 15 万个值作为参数展开，超出调用栈（`Maximum call stack size exceeded`）。改为逐个比较后 5.6 秒打开，并加测试防止模板再出现这种写法。
  - 组装软件新增 **Reference mapping** 组：只写了比对或一致序列软件（BWA、MIA、bcftools、iVar、GATK、Ion Torrent variant caller 等）的记录，多为人类和古 DNA 线粒体基因组（37,394 条中 35,226 条属哺乳纲）。补充 Sanger 编辑拼接软件（SeqScape、CodonCode、BioEdit、UGENE、ContigExpress、Phrap、DNAMAN、ChromasPro）、Minia、IVA、SMRT Analysis、MitoMaker、Supernova 和常见拼写错误（gnenious → Geneious、novaplasty → NOVOPlasty、mito-z → MitoZ）；MAFFT、MEGA、leeHom、HaploGrep 等不算组装。"Other" 从 18,353 条降到 1,699 条。不含字母的条目（"x"、"5.0"、"build 2025-03-24"）不再计为 Other。
  - `g2t-mitostats --exclude TAXON`（可重复）：去掉某个种（含亚种）或属、科、目、纲、门。*Homo sapiens* 占全部后生动物记录的 46%（68,433 条）。
  - 分类阶元缓存改放在数据目录（`gb/mitogenome-header/taxonomy_lineage.json`），换输出目录（例如加 `--exclude`）时不再重新联网查询。
  - 实测：首次下载 68 分钟、0 条失败；同日重跑 1 分 15 秒（150,025 条全部复用）。
- **分类核验评估**（`docs/evaluation/rdna-regions.md`）：在 Enchytraeidae（10,773 条）上检验按 feature 注释和按序列（Infernal cmscan + Rfam 5.8S/SSU/LSU 模型）核验基因归属。1,886 条 ITS 记录中 86% 的注释没有区段坐标；被同时归入 18s/28s 列的 1,485 条记录只含 18S/28S 引物端片段（最长 67 bp）。5.8S 在 1,883/1,886 条中定位成功。据此计划实现 `g2t-verify`。
- **修复：标题写作 "histone3 (H3)" 的记录未匹配**（Enchytraeidae 23 条）：补进 `gene_dict.yaml` 的 h3 同义词，并加测试。
- **文档**：README 改版（首页截图、功能表、快速开始、教程索引）；新增实例教程 `docs/tutorials/`（线粒体基因组统计、全后生动物、标本 × 基因矩阵），所有数字为实跑结果；`docs/CLI_REFERENCE.md` 加入两个新命令。

- **流水线失败时说明原因**：`StepResult` 新增 `error` 字段，extract、classify、organize 失败时填上原因；`g2t` 的失败信息从 `Step 1 failed` 变成 `Step 1 failed: no GenBank files (...) found in the inputs`，organize 写明缺哪几列（例：`Step 4 failed: missing required columns in <文件>: gene_type (<输出>)`）。`run()` 返回的 `log` 最后一行同此。
- **`g2t.extract`、`g2t.classify` 不再重新导出下划线开头的内部名字**（`_phrase_pattern`、`_yaml_path`、`_packaged_gene_dict`、`_load_gene_dict_yaml`、`_DEFAULT_*_GENES`、`_clean_outputs`、`_write_all`、`_process_record_loop`）；需要时从 `g2t._classify.<模块>`、`g2t._extract.<模块>` 导入。两个模块新增 `__all__`，列出文档里的接口（步骤函数、`main`、`MatchConfig`、词典函数、`TaxonomyService`、找 GenBank 文件的函数）；其余不带下划线的名字仍可导入，旧代码不受影响。测试检查两个模块不再带出内部名字。

- **拆分 `classify.py`**（1639 行 → 406 行）：`g2t/classify.py` 只留步骤函数 `classify()`、命令行和日志设置，其余代码原样移到 `g2t/_classify/`：`genedict`（同义词、用户词典、关键词正则）、`config`（`MatchConfig`）、`text`（短语匹配、DEFINITION 清洗、长度）、`tables`（记录过滤、输出文件名、列顺序）、`match`（第 1 轮）、`conflict`（第 1.5 轮冲突检查）、`recheck`（第 2 轮）、`report`（分基因文件、未匹配报告、record_status.csv）、`flows`（各轮在整张表上运行）。原来 `from g2t.classify import ...` 的名字照常可用，日志名仍是 `g2t.classify`。参数名、凭证清洗都没动。
- **修复：`python -m g2t.classify` 的 `classify.log` 是空的**（拆分时比对日志发现）。这种运行方式下模块名是 `__main__`，日志不进 `g2t` 的文件 handler。`classify.py`、`extract.py` 的日志名改为固定的 `g2t.classify`、`g2t.extract`；`g2t-classify` 命令本来就正常。

- **拆分 `extract.py`**（2137 行 → 383 行）：`g2t/extract.py` 只留步骤函数 `extract()` 和命令行 `main()`，其余代码原样移到 `g2t/_extract/`：`files`（找文件、打开、校验）、`parse`（一条记录 → 元数据行和 assembly 行、列顺序）、`taxonomy`（ete3 / Entrez 分类和缓存）、`status`（记录与文件状态、汇总、报告）、`process`（流式 / 内存读取、写表）、`batch`（每文件一进程、合并）。原来 `from g2t.extract import ...` 的名字都照常可用。日志名仍是 `g2t.extract`。示例数据上流水线 20 个 CSV、`--batch` 模式 10 个 CSV 与改动前逐字节相同，日志相同。
  - 唯一需要注意的：要替换分类服务（测试里的做法），改 `g2t._extract.taxonomy._taxonomy`；`g2t.extract._taxonomy` 不再导出，免得替换了却不起作用。

- **`g2t-extract` / `g2t-classify` 改为调用 `extract()` / `classify()`**（第四轮结构审查核对后）：原来命令行和函数各写一遍流程，修一个行为要改两处。现在命令行只解析参数，流程只在函数里。函数增加对应命令行选项的关键字参数（extract：表名、分隔符、`only_assembly`/`only_metadata`、`recursive`、`log_interval`、`update_taxonomy`、`report`/`details`/`tracking`；classify：`multiple_output`、`conflicted_output`、`moleculetype`、`mol_type`、`length`、`organelle`），默认值等于流水线原来的行为。输出不变：示例数据上两个命令和整条流水线的 CSV 与改动前逐字节相同。行为变化：
  - 失败时退出码为 1（原来 `g2t-extract` 所有文件都失败也返回 0；`g2t-classify` 报 Python 异常栈退出）。
  - 失败原因连同异常栈写进日志（原来流水线只剩 `Step N failed`）。
  - `extract(only_assembly=True)` / `only_metadata=True` 以各自的表为产物，判为成功（原来函数里没有这两种模式）。
  - `g2t-classify` 也写 `unmatched_report.txt`（原来只有函数写），每次运行先清掉旧的 `multiple_matched.csv`、`conflicted_genes.csv`（流水线同）。

- **版本名统一为 v0.04**：`g2t --version`（各子命令同）、`docs/CLI_REFERENCE.md` 标题输出 `v0.04`，`g2t_run.json` 新增 `"release": "v0.04"`；包版本仍为 `0.0.4`（pip/bioconda 用，PEP 440 会把 0.04 规范成 0.4），`g2t.RELEASE` 由它推出，两者不会再分开。
- **`--mtgenome_includes` 按给定的列执行**：原来正好给出旧默认的 6 列（coi 16s 12s cob cox2 cox3）时，会被当作“没给”，悄悄扩成全部 15 个线粒体列。默认行为不变：线粒体基因组复制到全部线粒体基因列，叶绿体基因组复制到全部叶绿体基因列（两份默认名单改为从 `g2t/genes.py` 生成）。README/USAGE 原来写成只复制到 6 列，已改正。
- **curate 在全部基因列里找登录号**：`DEFAULT_GENES` 改为 `GENE_ORDER`（原来只有 13 列，nd、atp 和叶绿体列中的登录号在校正表按登录号匹配时找不到）。
- `MatchConfig.gene_length_ranges` 删掉 17 个 `"none"` 条目（未列出的类型本来就是 none）；`ncbi_genes.MITO_GENES` 等写明只是预检计数用的常用标记。USAGE“使用现有环境”一句补上 pyyaml。

- **基因类型名单收成一份 `g2t/genes.py`**（第二轮结构审查核对后）：线粒体 / 核 / 叶绿体三组类型和矩阵列顺序（`MITO_TYPES`、`NUCLEAR_TYPES`、`PLASTID_TYPES`、`GENE_ORDER`）、凭证里的基因名（`VOUCHER_GENE_TOKENS`）。以下各处原来各写一份，现在都从这里导入：classify 的默认基因集合和 `--mtgenes_list/--ntgenes_list/--cpgenes_list` 默认串、organize 的 `gene_order` 和 `known_mito`（去掉了 yaml 里不存在的 `cytb`）、reconcile 和 curate 的 `GENE_TOKENS`。curate 原来那份缺 ND4L、EF-1 和叶绿体基因（`matK KUN 12345` 原来不剥掉 matK）。测试检查三组类型正好等于 `gene_dict.yaml` 的键，并且各步确实使用这份共享名单。
- **文档**：pyyaml 列入必需依赖（原来放在“可选依赖”下面却标着“必需”），“仅安装依赖”的 pip 命令补上 pyyaml；README 的 `--gene` 名单补 ND3、ND4L、ND6、ATP8；测试数量改为本次实测（589；3.12 与 3.9）。`docs/CLI_REFERENCE.md` 中 6 个命令的用法行写成了 `usage: make_cli_reference.py`，生成脚本已修好并重新生成；classify 默认基因串的顺序随 `GENE_ORDER` 变化（集合不变）。

- **基因词典只剩一份，分类器不再改模块全局**（据两份结构审查核对后修改）：
  - 删除 `classify.py` 里约 580 行的 `_BUILTIN_GENE_DICT`（与 `data/gene_dict.yaml` 逐条相同，41 类、同义词无差异）；pyyaml 改为必装依赖（这份内置副本原本就是为未装 pyyaml 准备的）；yaml 读不到时报错，不再静默换用另一份。
  - 删除 `GENE_DICT`、`_COMPILED_PATTERNS`、`_precompile_patterns`、`_reset_patterns`：词典放进 `MatchConfig.gene_dict`（每个配置一份），正则改为纯缓存 `_phrase_pattern`。`load_path_dict(path, base=None)` 返回合并后的新词典（原来就地修改全局，同一进程里后续运行都会带上）；新增 `default_gene_dict()`、`read_gene_dict_file()`。Python API 现在可以 `MatchConfig(gene_dict=load_path_dict("my.txt"))`。
  - `--path_dict` 文件不存在、格式不支持或没有条目时报错停止（原来只记一条 error，然后用默认词典照常分类）。`--mtgenes_list` 等列出的基因在加载 `--path_dict` 之后再检查（原来先检查，用户词典里的新基因会误报“不在词典中”）。
  - `g2t-classify` 每行日志只打印一次（原来 `_cli` 和 `setup_logging` 各挂一个控制台 handler，每行两遍）；完整日志写到输出目录的 `classify.log`，不再在当前目录写 `extract_methods.log`，也不再改根 logger 的级别。
- **三份基因名单对齐**：下载补上分类已有的 ND3、ND4L、ND6、ATP8（原来 `g2t-download --gene ND6` 报 Unknown gene），别名补 NAD2–NAD6、NAD4L；reconcile 的凭证基因前后缀补 ND4L、EF-1 和叶绿体基因（`KUN 12345 matK` 原来得不到凭证主体）。新测试 `tests/test_gene_dictionary.py` 检查三份名单一致，以后漏改一处测试即失败。RPB2、5.8S 只能下载、不单独分类，已在测试里列明。
- `reconcile` 不再引用 `voucher._is_blank_value`：空白判断移到 `utils.is_blank_value`（`voucher._is_blank_value` 保留为别名）。

- **读取 gzip 压缩的 GenBank 文件**：`.gb.gz`、`.gbf.gz`、`.gbff.gz` 等直接作为输入（文件发现、校验、计数、解析都支持）；示例数据压缩后得到的矩阵与未压缩时完全相同。
- **命令行默认流式提取**：`g2t` 与 `g2t-extract` 默认 `--stream`，与 Python API 一致（示例数据上两种模式的基因表与矩阵逐字节相同）；`--no-stream` 恢复内存模式，旧写法 `--stream` 照常可用。
- **矩阵新增 `type_material` 列**：模式标本信息（如 "holotype of …"）原本已随 source 限定词提取，但没有进入矩阵；Phyllodocida 有 570 条记录带此信息。
- **文档**：每个命令行参数都有说明（原 54 个没有，共 145 个），`docs/CLI_REFERENCE.md` 由 `scripts/make_cli_reference.py` 从各命令的 `--help` 生成；USAGE / USAGE.en 的输出文件说明列出全部中间文件及其内容。
- 代码注释去掉旧审查编号（P0-4、L-06 等）和"改进历史"式模块头，只保留对代码本身的说明。
- **规则冻结（frozen-heldout-v1）后、根据独立测试集和盲审发现的修正**：
  - 同一批次按 isolate 拆开的个体之间判为冲突，不再按“凭证号相同”合回（盲审 100 个合并中 5 个错误全属此类：MCZ 70173 个体 5/7；台纸 MW0003970 的 MW1-1/MW1-2）。
  - 编号主体两侧的数字段也比较：两边各有对方没有的编号 → 不同标本（`NSJC66_1005_S001` / `NSJC66_962_S001`）。
  - specimen_voucher 只有机构缩写、没有数字（`BRLU`）视为没有凭证号（isolate、单倍型标签不受影响）。
  - 间隔区记录同时写了基因时两列都放（“trnL gene…; and maturase K (matK) gene” → trnl-trnf,matk；“rbcL-atpB spacer…; atpB gene” → atpb,rbcl-atpb）。
  - 未注释叶绿体基因组 “UNVERIFIED: … chloroplast/plastid sequence”（≥20 kb）→ plastome；“gene for LSU/SSU rRNA”、“26S/25S rRNA”、“18S small ribosomal RNA” 写法补入；ITS 记录的 26S/25S 侧翼放入 28s 列。
- **审稿意见（R1–R3）标本连接修复**：
  - 编号主体相同，但两边都有不同的机构代码 / 采集人姓名（`USNM 123456` / `MCZ 123456`，`Christenhusz 5678 (K)` / `Hassemer 5678 (FLOR)`，`Juslen s.n. 2003` / `Muola s.n. 2003`）或不同的样品编号（`Polychaeta1D9_St95` / `Polychaeta1G8_St95`，同一站位的两个样品）→ 不同标本。一边没有这类词时照旧可合并（`WS399` / `COI_ZMMU_MSU_WS399`；`H & W 16149` / `C.B. Hellquist, J.H. Wiersema 16149`）。
  - 占位凭证号（`n/d`、`FCO:n./d.`、`s.n.`、`no voucher`、`NA`、`unvouchered` 等，且不含数字）视为没有凭证号，改用 isolate 等字段。
  - 一个凭证号是一个批次 / 一张标本台纸：同一凭证号下不同 isolate 测了同一基因（MCZ 70173 的 7 条 COI）→ 每个 isolate 一行，`match_evidence` 注明 “one voucher, several individuals”。
  - 矩阵物种名取该标本记录中出现最多的名字（原为字母序第一个）；名字不一致时新列 `organism_names` 列出全部。
- **下载默认排除转录组组装（TSA）**：`NOT "tsa"[Keyword]`，`--include-tsa` 加回；预检报告显示 TSA 条数（Phyllodocida 188,770 条中 162,319 条是 TSA）。
- `g2t_run.json` 记录每个输入 GenBank 文件的 SHA-256 和分类信息来源；openpyxl 改为必装依赖；不再提示缺少 pyyaml（内置词典与 yaml 相同）。
- **审稿意见（R2 M5）植物修复**：① mtgenome 同义词去掉过宽的 “complete sequence”“complete genome”（银杏 5 kb 的 rrn26+trnY 片段曾被当作完整线粒体基因组并复制进全部线粒体列），改为 “mitochondrion, complete sequence”“mitochondrial DNA, complete genome” 等；多基因记录只有标题写明 genome 才归为 mtgenome，否则保留各基因并标记冲突；② 第二轮 “circular” 先看细胞器：叶绿体 → plastome（银杏 AB684440 原被判为 mtgenome），并补 “chloroplast DNA, complete genome”“organelle: plastid”；③ 植物记录（Viridiplantae 或叶绿体组）第二轮丢弃其它细胞器的基因（叶绿体 16S/23S 不再进动物线粒体 16s 列；线粒体 rrn26 不进 28s）；动物保持原设计（“18S…; mitochondrial” 仍归 18s 并标记）；④ 新增蕨类常用 `rps4`、`atpb`。动物 5 个数据集归类结果不变（逐条对比）；Ophioglossaceae 新增 61 条 rps4、44 条 atpB。
- **再防两类错误合并**（Phyllodocida 复查）：① 尾号按原始分隔符判断，`R1941-END07-19-1` / `-2` 是批号（原来去掉分隔符后 “END07191” 被误当作第二个编号）；② 两组都有 specimen_voucher 且毫无共同编号时判为不同标本，即使 isolate 相同（`RUMF:ZO-00076` / `ZO-00070` 共用采集事件号 CF21112021）。
- **批号规则细化**：一边无尾号、另一边尾号含字母（`28S_ZMMU_MSU_WS16799_Pr21` / `COI_ZMMU_MSU_WS16799`，Pr21 为附加样品号）不再视为不同标本；只有两边尾号不同（`-E01`/`-E02`）或一边多出纯数字尾号（`CE438`/`CE438_2`）才算。Priapulidae 4 对因此恢复合并。
- **植物（叶绿体）基因**：按动物的逻辑增加第三组“叶绿体”：`/organelle` 为 plastid/chloroplast，或标题含 chloroplast/plastid（不含 “nuclear gene for chloroplast product”）。新增 17 类：`plastome`（叶绿体基因组，优先匹配；>2 个基因且 ≥20 kb 也归为 plastome，登录号复制到各叶绿体基因列，同 mtgenome）、`rbcl`、`matk`、`trnl-trnf`、`psba-trnh`、`rpl16`、`rps16`、`ndhf`、`rpoc1`、`ycf1` 及间隔区 `rpl32-trnl`、`trns-trng`、`atpf-atph`、`petb-petd`、`psbe-petl`、`rbcl-atpb`、`trnk`（内含子）。间隔区优先于基因（“rbcL-atpB intergenic spacer” 不算 rbcL 基因；基因同义词只用 “(rbcL)”“rbcL gene”、全称）。叶绿体记录按自己的长度范围 150–300,000 bp（`--length2_cpgenes`），叶绿体基因组不再被全局 50 kb 上限过滤；动物数据结果不变。植物核 26S/25S rRNA 归入 `28s`；植物线粒体写法 NADH1–NADH6 归入 nd1–nd6。`g2t-download --gene` 支持 RBCL、MATK、TRNL-TRNF、PSBA-TRNH 等。实测（2026-10-23 下载）：Ophioglossaceae 归类记录 143 → 4,067、≥2 个基因的标本 3 → 1,023；Nymphaeaceae 772 → 1,932、352 → 644；Ginkgoaceae 82 → 325、27 → 130。
- **修复**：reconcile 合并后的新名字只取最高优先级编号字段的编号主体，该字段无编号（如 voucher 写 “s.n.”、编号只在 isolate）时崩溃（Nymphaeaceae 上发现）。
- **新增线粒体基因类型 nd1–nd6、nd4l、atp6、atp8**（共 22 类）：原来这些记录全部未归类（Hesionidae 284 条未归类记录中 221 条是 ATP6/ND4/ND5）。同义词含 “NADH dehydrogenase subunits 4”、“ATP synthase Fo subunit 6” 等实际写法；核基因 “ATP synthase subunit beta” 不归入。矩阵新增 9 列（cox3 之后），线粒体基因组的登录号也复制到这些列。
- **下载断网不再退出**：NCBI 请求失败最多重试 10 次（等待 2、4、8… 秒，单次不超过 120 秒，合计约 8 分钟），原来 5 次共约 30 秒，Phyllodocida 下载中一次 `RemoteDisconnected` 即退出。
- **同一凭证号、物种名写法不同 → 同一标本（一行）**：如 SIO:BIC:A10920 的 18S/28S/线粒体基因组登记为 *Bathykurila* sp.，H3 登记为 *Bathykurila* sp. AH-2023，原来分成两行。现在凭证号相同（忽略大小写和标点，如 `SIO:BIC:A6302` = `SIO-BIC A6302`）且属名相同即合并，`match_evidence` 注明 “organism names differ”。属名不同（WS3020：*Halicryptus* / *Priapulus*）仍为两行，两行的 `match_evidence` 写 “same voucher as … (different genus, not merged; check)”，并打印警告。
- **两组都有同一基因时，“同一篇论文”不再足以合并**：CE438 与 CE438_2（Martin et al. 2024 列为两个标本，各有一条 28S）原来因同一篇论文被并成一行（28S 格两个登录号）。现在需要同一采集日期或同一坐标（同一凭证号在这种情况下也只算中等证据）；未合并的两组在 `match_evidence` 写 “similar voucher to … (both have the same gene, not merged; check)”。
- 文献表基准（`benchmark/`）：Polynoidae 33 个标本，标本连接召回 86.6% → 96.2%（精确度 100%）；Haplotaxidae 75 个标本，召回 88.1% → 90.1%。
- **跨区段 rRNA 记录放进它覆盖的每一列**：标题同时含 ITS 和 18S/28S 的记录（如 OM007981 “small subunit ribosomal RNA gene and internal transcribed spacer 1”，1908 bp；OM105844 “internal transcribed spacer 2 and large subunit ribosomal RNA gene”）原来只进 `its1-its2` 列，现在 `gene_type` 为 `18s,its1-its2` / `28s,its1-its2`，矩阵中同一登录号同时出现在 18s（或 28s）列和 its1-its2 列，`Assignment_reason` 注明 “multi-region rDNA”。完整 18S–28S 记录仍为 `18-28s`（organize 已复制到 18s/28s/its1-its2 列）。Polynoidae 文献表中 60 个 18S/28S 登录号因此出现在矩阵的 18s/28s 列。
- **输入扩展名**：`g2t`/`g2t-extract` 现在识别 `.gbf`（NCBI 及很多工具导出的 GenBank 文件扩展名）；全部支持 `.gb`、`.gbf`、`.gbk`、`.gbff`、`.genbank`（不区分大小写）。此前 `.gbf` 文件被跳过，报 “no GenBank files found”。

## v0.04（2026-10-22）

包内版本 `0.0.4`。本版重写下载步骤（按记录增量、`--since`、`--no-store`），并修复全库代码审查发现的问题。

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
- 命令行整理：选项按 selection / update / run 分组，帮助文字简短；`--store` 与 `--no-store` 互斥；`--since`、`-w`、`-b`、`--minlen/--maxlen` 在解析时校验（错误退出码 2）；新增 `-n`（= `--dry-run`）、`-q`（不显示进度）；进度输出到 stderr；`--resume` 隐藏（续传自动）。
- 测试：新增 24 项，共 412 项通过。

### 修复（2026-10-22 全库代码审查）
- **`g2t` 成功时退出码为 1**：`main_pipeline` 返回了 `PipelineResult` 对象（`sys.exit(obj)` → 退出码 1 并把对象打印到 stderr）。现返回 0。BUGS.md 中 2026-05-12 的条目原因已更正。
- **下载**：`--gene 12S/16S` 的 Entrez 子句缺少括号（`A OR B AND mitochondrion[filter]` 被解析为 `(A OR B) AND ...`，漏掉无 `/organelle` 的记录；Annelida 16S 少 276 条）；`--since` 只补缺失的 accession.version，记录内容更新但版本号未变（如补注文献、改正产地）时不重新下载，现重新下载；选择变为空时目录不刷新（旧批文件和 manifest 保留）；“未变化”快捷路径只数批文件个数不看文件名，`--no-store` 留下编号空洞后会误删最后一个批文件；分页用历史检索自身的 Count。
- **分类**：`18-28s` 的 5 条同义词写成 “5.8S rRNA gene”，而 NCBI 标题是 “5.8S ribosomal RNA gene”，整段 rDNA 记录全部落到 `its1-its2`；核糖体冲突检查把 `12s,16s`/`18s,28s` 双匹配截成单个基因；只看 `/organelle` 分组，无 `/organelle` 但标题含 mitochondrial 的 16S 记录按核基因处理成 `28s`；`mtgenome` 长度区间设为 `none` 时仍强制 3000 bp，`_conflict_mito` 写死 3000；`.py` 词典键未小写；16S 甲基转移酶排除式不认 `16S rRNA (guanine(527)-N(7))-methyltransferase`。
- **提取**：流式模式的列预扫描以系统编码打开文件，遇到一个 Latin-1 字节即静默放弃，lat_lon/collection_date/Ref* 等所有非基础列丢失而状态仍为 SUCCESS；同一记录出现在两个文件时装配表左连接使行数翻倍；`--log_interval 0` 除零；批处理在 spawn/forkserver 启动方式下子进程无分类学后端（Class/Order/Family/Genus 全空）；注释自由文本中同一行多个装配字段互相吞并；`os.walk` 子目录顺序随文件系统变化、非递归模式把子目录当文件。
- **凭证**：`reconcile` 合并边的顺序随输入行序变化，冲突时得到不同的聚类；纯标点凭证（`..`、`-`）未像 `voucher` 那样回退到 isolate，导致永不参与合并；`organize` 的 `--metadata_columns` 含空格/括号/关键字的列名经 `itertuples` 改名后取不到值（输出空）；`curate` 中矩阵里很短的凭证（如 `12`）会匹配任何包含该 token 的用户凭证，校正写到错误的行；登录号不在矩阵中时错误信息误导；`curated_fields` 为 NaN 时写成 `nan; lat_lon`；凭证列（`USNM123456`）被当成基因列。
- **pipeline**：`_find_fallback` 找到同名候选文件后什么都不做，下一步仍读固定路径失败；现直接报错并提示候选文件名。`g2t` 删除从未生效的 `--*_extra` 参数。
- **scripts**：`download_edirect.sh` 用了 EDirect 不存在的 `-retstart/-retmax`（每批都失败），改用 `-start/-stop`；`generate_ncbi_update.py` 的 geo_loc_name 列键名错误（输出表永远为空）、空坐标写成 `nan S nan W`、空国家写成 `nan`、全称月份未缩写、日期区间丢掉结束日、单个词的凭证抛 IndexError、final.csv 无 isolation_source 列时对所有记录填充；`validate_against_excel.py` 准确率分母用错（出现 346/345 = 100.3%）、空凭证按 `"nan"` 比较；`download_by_accession.py` 失败批次不记录，`--resume` 永不重下；`download_requests.py`/`download_by_accession.py` 不校验批内记录数，截断响应被当成功；`download_requests.py` 不发送 email；`compare_phyllodocida.py` 在 Python 3.9 导入失败（运行时 `str | None`），并删除写死在源码里的 NCBI API key（改读 `NCBI_API_KEY`）；`build_gene_dict_v2.py` 不检查 SynGenes mt 目录。
- 测试：新增 24 项回归测试，共 436 项通过。
- 多版本验证：在 Python 3.9.23（pandas 2.3 / biopython 1.85）、3.10、3.12、3.13.16（pandas 3.0 / biopython 1.88）下 436 项测试全部通过，ruff 与 mypy 均通过。去掉 `infer_objects(copy=False)`（pandas <2.0 无此参数，3.0 起弃用）；dev 依赖加入 `types-PyYAML`（无它时 mypy 在 3.9 下报 yaml 缺 stub）。
- 代码整理：`ruff check src/ tests/`（CI 的 lint 步骤）与 `mypy src/` 现在全部通过（原 354 条 lint、41 条类型错误）：类型注解统一为 PEP 585/604 写法（所有模块已有 `from __future__ import annotations`，Python 3.9 可用）、删除未用导入/变量、隐式 Optional 补全、Biopython 注解值统一 `str()`。删除从未可达的装配提取 “strategy 3”（`extract_assembly_from_raw_block`）和 `compare_results.py` 未调用的 `parse_docx`；`merge_literature_metadata.py` 经纬度改为 GenBank 的 `57.68 N 11.96 E` 格式；两脚本重复 accession 统一为首行优先；`format_voucher` 的 “RBINS IG 33953” 分支改为 `RBINS:IG:33953`（institution:collection:id）并置于通用分支之前。

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
