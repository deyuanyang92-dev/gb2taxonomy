# 教程：标本 × 基因矩阵（download → g2t → curate）

[← README](../../README.md) · [线粒体基因组统计](mitogenome-stats.md) · [全后生动物](metazoa.md)

目标：把一个类群在 GenBank 里的 COI、16S、18S、28S 等记录，按标本（物种 + 凭证号）合成一行一个标本的矩阵，供多基因系统发育使用。

下面第 1 步要联网；第 2、3 步用仓库自带的离线数据 [`examples/haplotaxidae`](../../examples/haplotaxidae)（347 条记录，来自 Martin et al. 2024 *Zool J Linn Soc* 表 S1），任何电脑上都能复现，输出数字是 2026-11-13 实跑结果。

## 1. 下载（联网）

先看 NCBI 里有什么，不写文件：

```bash
g2t-download -t Priapulidae -o gb --dry-run
```

报告该类群的记录构成（WGS、mRNA、TSA……）、各基因记录数、要下载多少条。确认后正式下载：

```bash
g2t-download -t Priapulidae -o gb             # 默认 "markers"：不含 WGS/mRNA/TSA/RefSeq
g2t-download -t Priapulidae -o gb --gene COI,18S,28S   # 只要某几个基因
g2t-download -t Priapulidae -o gb --mito                # 只要线粒体记录
```

2026 年 10 月实测：Priapulidae（txid37891）共 95,726 条 nuccore 记录，其中 69,704 条 WGS contig、23,053 条 mRNA；默认选择留下 548 条标记基因记录。

输出目录：

```
gb/
├── _records.sqlite          所有下载过的记录（按 accession.version 存）；重跑只下新的
└── markers/
    ├── batch_0001.gb ...    GenBank 文件
    ├── accessions.tsv       本次选择的登录号
    ├── changes.tsv          与上次相比：新增 / 更新 / 撤回
    └── manifest.json        查询式、日期、完成状态
```

以后更新：同一条命令再跑一次即可，只下载新增或版本号变了的记录（`--since auto` 只查上次以来修改过的记录，更快，但看不到被撤回的）。

## 2. 生成矩阵（离线示例）

```bash
g2t -i examples/haplotaxidae/haplotaxidae.gb -o hap --stream
```

约 4 秒。结果：347 条记录全部分到基因类型，矩阵 `hap/organized_genes/organized_species_voucher.csv` 共 95 个标本，其中 70 个有 ≥ 2 个基因。例如标本 CE1233（*Haplotaxis* sp. SE）一行：

| organism | coi | 16s | 12s | 18s | 28s |
|---|---|---|---|---|---|
| *Haplotaxis* sp. SE | PP988405.1 | PP994580.1 | PP968574.1 | PP990634.1 | PP994523.1 |

各步的中间文件和用途见 [README › Steps 1–4](../../README.md#steps-14--the-pipeline-g2t)。要核对某条记录为什么没进矩阵，看 `hap/labeled_genes/record_status.csv`（每条记录一个状态和原因）。

## 3. 用自己的表校正（g2t-curate）

GenBank 里常缺坐标、采集日期，或物种名已过时。把你自己的表（任何列名，中英文都行）合进矩阵：

`my_table.csv`（**示例数值，仅演示用法，不是真实采集记录**）：

```
登录号,纬度,经度,采集日期,备注
PP988405,59.85,17.63,2019-06-12,补自原文表 S1
PP988416,45.2,9.1,,
```

```bash
g2t-curate -m hap/organized_genes/organized_species_voucher.csv -u my_table.csv -o hap/curated.xlsx
```

```
Curated: hap/curated.xlsx (4 changes in 2 specimens; 0 problems)
  登录号      -> accession
  纬度        -> latitude
  经度        -> longitude
  采集日期    -> collection_date
  备注        -> (added as a new column)
```

`curated.xlsx` 的 *Changes* 表逐条记录改动（旧值 → 新值、按什么匹配的）：

| specimen | field | old_value | new_value | matched_by |
|---|---|---|---|---|
| Haplotaxis_sp_SE_CE1233 | lat_lon | | 59.85 N 17.63 E | accession |
| Haplotaxis_sp_SE_CE1233 | collection_date | | 12-Jun-2019 | accession |
| Haplotaxis_sp_SE_CE1233 | user:备注 | | 补自原文表 S1 | accession |
| Acanthobdella_peledina_CE27296 | lat_lon | | 45.2 N 9.1 E | accession |

列名认错了就用 `--map "编号=voucher"` 指定。不想自己建表，可先生成预填模板，改格子后再合并：

```bash
g2t-curate -m hap/organized_genes/organized_species_voucher.csv --template corrections.xlsx
g2t-curate -m hap/organized_genes/organized_species_voucher.csv -u corrections.xlsx -o hap/curated.xlsx
```

期刊信息过时（GenBank 仍写 "Unpublished"）可以用 `g2t-cite` 自动找，结果表的 *Curate* 页就是 `-u` 能直接用的表，见 [线粒体基因组教程 › 补期刊](mitogenome-stats.md#4-补上已发表但-genbank-没更新的期刊g2t-cite)。
