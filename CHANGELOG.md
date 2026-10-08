# 更新日志 / Changelog

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
