# 评估：污染与鉴定错误提示（`g2t-screen`，v0.05）

[← rRNA 区段评估](rdna-regions.md) · [← README](../../README.md)

## 问题

GenBank 里有些序列的物种名是错的：样品污染（混入别的动物的 DNA）、标本鉴定错误或样品混淆。建树时它们表现为某个种落进另一个属里。

过去的检查做法是：比对 → 算遗传距离 → 建树 → 综合判断，全靠手工，费时。`g2t-screen` 把“算距离、找异常”这一步自动化，对每个标记（COI、16S、18S……）分别检查，几秒到一分钟出结果。

## 方法

不需要外部程序，不做多序列比对。

1. 每条序列取正反两条链的全部 k-mer。这样，反向提交的记录和含内含子的 COI 也能正确比较。
2. 用 5-mer 组成的余弦相似度为每条序列找出 20 个候选近邻。
3. 用 11-mer 的 Mash 距离 −ln(2J/(1+J))/k 计算到候选近邻的距离。这个距离在 ≤ 5% 时与 p 距离很接近。
4. 给出两类提示：

| 提示 | 条件 | 适用标记 |
|---|---|---|
| `near_other_genus` | 最近的序列属于另一个属，且距离 ≤ 2%；同属最近的序列至少再远 2 个百分点，或同属没有别的序列 | 进化快的标记：COI、16S、12S、CYTB、ND 系列、ITS 等。18S、28S、H3 在不同属之间可能完全相同，不按属判断 |
| `same_species_far` | 同名（已定名种）的两条记录相差 > 15%；注释同时给出各自到本属的距离，离本属远的那条更可疑 | 快速进化的标记，含线粒体基因组中的基因 |
| `remote` | 本数据集中同一标记里，没有任何序列的 8-mer 距离在 0.35 以内（约相差 30% 以上） | 所有标记，序列数 ≥ 10 时才判断 |

提示只是“值得核查”，不等于结论。可能的原因还包括：属本身不是单系、名称过时、近缘属之间共享单倍型。

## 结果（2026-11-13）

| 数据 | 检查的序列数 | 被提示 | 用时 |
|---|---|---|---|
| Opheliidae COI | 293 | 1 | 5 秒 |
| Hesionidae 全部标记 | 2,071 | 0 | 19 秒 |
| Enchytraeidae 全部标记 | 9,188 | 10 | 57 秒 |

**Opheliidae**：只提示 1 条，就是此前人工建树发现的那条。

| 登录号 | 名称 | 提示 |
|---|---|---|
| OQ417162.1 | *Polyophthalmus pictus* | 与一条 *Armandia* 序列相差 0.4%；最近的 *Polyophthalmus* 相差 27.9% |

**Enchytraeidae**：提示 10 条，均待原文核对。

| 登录号 | 标记 | 名称 | 提示 |
|---|---|---|---|
| ON419122.1 / ON419111.1 | COI | *Fridericia loretensis* / *Hemienchytraeus tenuiculus* | 两条 COI 相同（0.0%），互为最近；与各自同属最近的序列分别相差 18.1% 和 16.8%。两条中至少有一条名称有误 |
| GU901794.1 | 16S | *Enchytronia parva* | 与 *Enchytraeus* 相同（0.0%） |
| GU902056.1、MG252150.1 | COI | *Enchytronia parva* | 与 *Enchytraeus* 相差 1.2–1.5% |
| PP169244.1 | ITS | "*Enchytraeus parva* complex sp. L3" | 与 *Enchytronia* 相差 0.1%，与前几条反向对应 |
| PZ404487.1 | COI | *Enchytraeus* sp. | 与 *Achaeta* 相同 |
| PZ404491.1 | COI | *Marionina* sp. | 与 *Enchytraeus* 相差 1.1% |
| LN999380.1、LN999152.1 | COI | *Marionina* sp. 1 RV-2016 | 与 *Globulidrilus* 相差 0–0.4% |

在 `g2t` 流程中，第 5 步在第 2b 步之后运行，结果是 12 条。原因是：原来标为 "28s,its1-its2" 的记录改为纯 ITS 后，ITS 组的序列变多，多出 3 条 ITS 提示：

| 登录号 | 名称 | 提示 |
|---|---|---|
| PP169655.1 | *Marionina argentea* complex sp. CED | 与 *Cernosvitoviella* 相差 1.4% |
| PP168965.1 | *Cernosvitoviella atrata* complex sp. CEB | 与 *Marionina* 相差 1.4% |
| PP955315.1 | *Globulidrilus riparius* | 与 *Marionina* 相差 0.4% |

另外，PP169244.1 不再出现在提示中。

几条“*Enchytronia parva*”记录与 *Enchytraeus* 相差 0–1.5%，而 ITS 记录又写作 "*Enchytraeus parva* complex"。这看起来是名称或归属问题，而不是污染（推测，未核对原文）。

## 更正：之前的原型把两条序列误判为“远离所有序列”

原型用 MAFFT 比对后算 p 距离，得出：

- AY364856.1（*Ophelia rathkei*）与所有序列相差 ≥ 48.7%；
- LC670764.1（*Travisia*）相差 51%。

用 k-mer 方法重算（正反链都考虑，不依赖比对）：

- AY364856.1 与另一条 *Ophelia* 相差约 25%；
- LC670764.1 与另一条 *Travisia* 相差约 21%；
- 两者都在科内属间的正常范围，因此不提示。

那两个约 50% 的距离，很可能是 MAFFT 处理反向或含内含子的序列时比对错位造成的。这也是不再用多序列比对来筛查的原因。

## 新增：同名记录差异过大（`same_species_far`，2026-11-13）

起因是 *Ophelia* 的 16S。研究中已经人工剔除了 HM746718（"*O. neglecta*"）和 PX620130（"*O. translucens*"），但原来的两类提示都没有发现它们。

**规则：** 两条记录的种名相同（已定名的种，不含 sp.、cf.、aff.、complex），在快速进化的标记上相差 > 15%，就提示这两条。

**配套改动：**

1. **距离改为按较短序列计算共享 k-mer：** 1 − C^(1/k)，C = 共享 k-mer 数 / 较短序列的 k-mer 数。原来的 Jaccard / Mash 距离会把片段与完整基因的差异算大：544 bp 的片段与包含它的 1,231 bp 16S 完全相同，Jaccard 也会算出 4.5% 的差异。
2. **线粒体基因组也参与比较：**
   - 蛋白基因从注释中提取；
   - 12S/16S 没有注释时，用 `g2t-verify` 的 HMM 定位。例如 OZ260669 完全没有注释，原来无法比较。

**Opheliidae（697 条标记记录 + 13 条线粒体基因组）**

| 登录号 | 名称 | 提示 |
|---|---|---|
| **HM746718.1** | *Ophelia neglecta* 16S | 与 *O. neglecta* 线粒体基因组 OZ260669 的 16S 相差 **24.3%**；与其他 *Ophelia* 都相差 ≥ 23%；最近的是 *Ophelina*（16.2%）。OZ260669 与 *Ophelia*/*Euzonus* 相差 7.8–8.5%，正常。**问题在 HM746718。** |
| **DQ779637.1** | *Travisia forbesii* 16S | 与 *T. forbesii* 线粒体基因组 OZ479805 相差 23.6%；最近的是 *Ophelina*（18.3%）。OZ479805 与 *Travisia* 相差 7.5%，正常。情况同上 |
| ON716040.1 | *Polyophthalmus pictus* COI | 与另一条 *P. pictus* 相差 18.8%（可能是隐存种，推测） |
| MN539268.1 | *Thoracophelia dillonensis* COI | 与 KC164682 相差 18.4%；与另一条 *Thoracophelia* 只差 0.2% |
| OQ417162.1 | *Polyophthalmus pictus* COI | 与 *Armandia* 相差 0.2%（`near_other_genus`） |

HM746718 和 DQ779637 都是较早的 16S 片段，都最接近 *Ophelina*，而不是自己所在的属（推测：可能是同一批样品或同一实验室的污染，未核对）。

PX620130 *O. translucens* 没有同名的其他记录，这条规则判断不了。它与所有 *Ophelia* 相差 21–23%，与 *Ophelina acuminata* 相差 14.6%，值得核查；要判断需要建树。

**其他科的提示量：** Hesionidae 17 条，Enchytraeidae 34 条。多数集中在 *Fridericia*、*Marionina argentea*、*Buchholzia fallax*、*Enchytraeus buchholzi* 等，推测多为隐存种复合体，即同一个名称下有多个种。这类提示本身也有用：说明这个名称的序列不能直接当作一个种使用。

## NCBI 自己存疑的记录（`ncbi_unverified`，v0.06）

GenBank 对无法核实的记录有两种标注：

- 标题以 "UNVERIFIED:" 开头；
- COMMENT 写 "GenBank staff is unable to verify source organism provided by the submitter"（来源物种无法核实），或 "… unable to verify sequence and/or annotation …"（序列或注释无法核实，例如 COI 翻译不出正常蛋白）。

`g2t-screen` 把这些记录单独列为 `ncbi_unverified`。

**结果：**

- Hesionidae 有 28 条。
- 按距离提示的 22 条中，有 6 条同时被 NCBI 标为存疑：
  - PX640331、PX640352、PX640417、PX640477（"UNVERIFIED:" COI，都与同名记录相差 > 18%）；
  - MW278311（来源物种无法核实；与同名的 *Leocrates chinensis* 相差 20.6%，最近的是 *Sirsoe*）；
  - PP652214。

两种完全独立的依据（距离、NCBI 人工审核）指向同一批记录，说明距离提示是有效的。

## 测试过但没有采用：“比本属更接近另一个属”

**规则：** 另一个属在 10% 以内，且比本属最近的序列至少近 3 个百分点，就提示。用全部 Opheliidae 标记数据（697 条）测试。

**结论：** 这条规则在 Opheliidae 提示 14 条、Hesionidae 11 条、Enchytraeidae 28 条，大多是 28S。在 28S 上，k-mer 距离无法按属判断：

- 28S 保守，属间常只差 2–4%；
- 各记录测的是 28S 的不同区段。有时“本属最近的序列”显示 100%，就是因为没有重叠区段。

Opheliidae 被提示的 28S 涉及 *Armandia*、*Travisia*、*Polyophthalmus*、*Ophelina* 等，不可能都有问题。所以这条规则没有采用。

**Ophelia 的情况（2026-11-13，697 条 Opheliidae 标记记录）**

- 111 条 *Ophelia* COI 都与本属序列最近。
- *Ophelia rathkei* 的三条记录都偏离本属：
  - COI：AY364856 离最近的 *Ophelia* 24.9%，离其他属 43%；
  - 28S：AY366513 离 *Ophelia* 13.5%，离 *Ophelina acuminata* 6.1%；
  - 28S：X80651 离 *Ophelia* 9.1%，离 *Thoracophelia* 3.6%。
- 这可能是鉴定问题，也可能是该种的归属问题（推测）。
- 要判断，需要在重叠区段比对后建树，即后续计划的“系统发育位置检查”。k-mer 距离做不到。

## 不做的检查

- **终止密码子。** *Travisia* 的 COX1 含 II 类内含子（Sci Rep 2022），按终止密码子判断会误报。
- **按属聚类、看大距离。** 属间距离超过 20% 时，距离已经饱和（*Armandia* 的原型误报就是这样来的）。所以只在距离很小时比较不同的属。
