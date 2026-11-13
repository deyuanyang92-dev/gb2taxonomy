# 评估：用基因树检查污染和鉴定错误（`g2t-screen --tree`）

[← 距离初筛评估](screen.md) · [← README](../../README.md)

## 为什么要建树

距离初筛（k-mer）很快，但实测有三处做不到：

1. **找不到外部来源。** PX620130 "*Ophelia translucens*" 的 16S 与 *Hyboscolex*（鳞腹沙蚕科 Scalibregmatidae）相差只有 0.2%，但海蛹科数据里没有 *Hyboscolex*，距离法无从比对。
2. **距离饱和。** 序列间相差 20% 以上时，“离本属远”和“离别属近”分不开。
3. **保守基因没法按属判断。** 28S 在属间常只差 0.2–4%。

做法参考 OrthoFinder 3（Emms et al. 2026，*Nat. Methods* 23:1327）：先宽收，再由系统发育决定归属。默认软件与它相同：FAMSA 比对，去掉空位 ≥ 90% 的列，FastTree 建树。

## 流程

1. **背景类群。** 沿 NCBI 分类往上找，直到找到一个包含至少 6 个其他科的祖先，再取其中标记记录最多的 6 个科。例如海蛹科（Opheliidae）→ 小头虫目一类（Scolecida）：Capitellidae、Orbiniidae、Maldanidae、Paraonidae、Arenicolidae、Scalibregmatidae。也可用 `--background` 指定。下载是增量的。
2. **取样。** 研究类群全部取；某个标记超过 800 条时，每个种取一条，再加上距离初筛提示过的序列。背景类群每个科最多 12 条，一个种取一条。线粒体基因组中的基因也加入。
3. **统一方向。** 每条序列与最长的几条（不含它自己）比较 11-mer，决定是否反向互补。FAMSA 不会自动调方向。
4. **比对与建树。** FAMSA（没有则用 MAFFT）→ 去掉多空位的列 → FastTree GTR+Γ（超过 1,000 条用 `-fastest`）。`--iqtree` 改用 IQ-TREE（ModelFinder + 1,000 次 UFBoot）。
5. **无根判断。** 对每条研究类群的序列，找出满足以下条件的最小一侧，作为它的“近邻”：
   - 来自一个支持率 ≥ 0.9 的分割；
   - 包含这条序列，另含至少 2 条覆盖良好的序列。

   这样不需要定根，结果不受“以哪个科作外群”影响。
6. **区分记录错误和分类问题。**
   - **记录错误：** 近邻全是别的科（一到两个科），或全是本科的别的属；而本属的其他序列不在这里。
   - **分类问题，写入 `phylo_taxon_notes.csv`：**
     - 本属的序列也在这里；或本属超过一半的序列都被提示；
     - 两个属互相落在对方里（异名，或属的界限未定，例如 *Ophiodromus*/*Oxydromus*）；
     - 这条序列是本属唯一的一条（只有当它与别科序列相差 ≤ 5% 时，才按混入提示）。
7. **不作判断的情况：**
   - 覆盖不到核心区段（大多数序列都有数据的那些列）的一半，说明测的是基因的另一段；
   - 近邻混有本科和其他科；
   - 18S、28S 不做属一级的判断。
8. **距离复核。** 只有当这条序列按 k-mer 距离也更接近它所在的那一组、而不是本属最近的序列时，才保留提示。高支持率的位置也可能是假象，见下文的 HQ025027。
9. **跨基因印证（按标本）。**
   - 同一标本的其他基因都在正常位置 → 这条序列本身有问题（污染或混样）；
   - 其他基因也错位 → 可能是标本鉴定错误。

## 结果（2026-11-13）

**Opheliidae**（6 个标记；297 条 COI、153 条 16S 等，另加背景序列；含下载约 3 分钟，只建树约 45 秒）

| 登录号 | 名称 | 位置 | 支持率 | 与所在组的距离 / 与本属的距离 | 判断 |
|---|---|---|---|---|---|
| **PX620130.1** | *Ophelia translucens* 16S | Scalibregmatidae（*Hyboscolex*） | 0.93–0.99 | 0.2% / 17.9% | 鳞腹沙蚕科序列，混入或鉴定错误 |
| **HM746718.1** | *Ophelia neglecta* 16S | Capitellidae（*Notomastus*） | 0.99 | 11.8% / 17.9% | 小头虫科序列；同一批提交（HM7467xx）中也有 *Notomastus*（推测为混样） |
| **OQ417162.1** | *Polyophthalmus pictus* COI | *Armandia* | 1.00 | 0.25% / 22.4% | *Armandia* 序列 |

这三条都在人工核查结果中：前两条是研究中从 MEGA 矩阵剔除的序列，OQ417162 是人工建树时发现的。

被距离复核排除的 4 条（树上的位置是假象）：

| 登录号 | 树上的位置 | 实际距离 | 原因 |
|---|---|---|---|
| HQ025027 *Travisia forbesii* COI | Cirratulidae（1.00） | 与另一条 *T. forbesii* 相差 3.9%，与 Cirratulidae ≥ 31% | 长枝吸引 |
| LC677172 *T. sanrikuensis* COI | Scalibregmatidae（0.95） | 与 *Travisia* 相差 12.5%，与所在组 22.6% | *Travisia* 本身的位置不稳 |
| HQ025072 *Ophelina acuminata* 28S | *Ophelia*（0.98） | 与另一条 *O. acuminata* 完全相同 | 该区段太保守，分不开属 |
| AY040681 *Armandia maculata* 18S | Capitellidae/Orbiniidae（0.94） | 与本属 6.4%，与所在组 7.9% | 位置不可靠 |

*O. rathkei* 的三条序列都没有被提示：X80651 与 *O. neglecta*、*O. bicornis* 聚在一起（0.93），另两条位置不确定。此前按距离对它的怀疑不成立。

**Hesionidae**（背景：Nereididae、Syllidae、Pilargidae、Chrysopetalidae、Phyllodocidae、Glyceridae；约 70 秒）

- **记录错误 2 条：**
  - DQ442580 *Podarkeopsis arenicolus* 16S，落在 Pilargidae（0.99）；
  - DQ442575 *Leocrates chinensis* 16S，落在 *Dalhousiella*/*Leocratides*（0.93）。
  - 距离初筛也曾把这两个名称标为同名记录差异过大。
- **分类问题 14 条，例如：**
  - *Ophiodromus* 与 *Oxydromus* 互相落在对方里（通常视为异名）；
  - *Paraleocrates* 落在 *Leocrates* 里（0.4%）；
  - *Dalhousiella* sp. T72992-1 的 COI、16S、28S 都落在 *Leocratides* 里，跨基因一致，因此是标本或属的问题，而不是单条序列的问题。

## 局限

- **需要外部软件：** FAMSA 或 MAFFT，加上 FastTree 或 IQ-TREE。安装命令：`conda install -c bioconda famsa fasttree`。
- **背景类群要包含真正的来源**，才能指出序列来自哪个科。来源不在背景里时，这条序列通常位置不确定，不会被提示，只能靠距离初筛中的“与所有序列都很远”（`remote`）。
- **FastTree 的支持率（SH-like）通常偏高。** 因此另加了距离复核，也可以用 `--iqtree` 取 bootstrap。
- **阈值只用了两个科来校准：** 支持率 0.9、覆盖率 0.5、离别科序列 5%。
