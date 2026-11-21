# 方案: 未写明组装/测序方法的记录 → 追踪原始文献 (2026-11-22, 待新会话执行)

## 已有基础
- `g2t-cite` (src/g2t/cite.py): 对 "Unpublished" 记录按 (标题, 作者) 去重, 查 Crossref, 打分 high/medium/check/low, 输出 xlsx, 缓存 crossref_cache.json。
- `mitostats`: assembler_group = "No assembly info in record" / "Not stated" / "Unknown"(平台) 即待追踪对象。

## 探查结果 (EBP 线粒体 6309 条)
- 1803 条 "No assembly info": 只有 7 条有文献标题 (3 个不同标题)。其余是基因组项目批量记录, 没有论文 → 追文献对这类无意义, 只能靠 COMMENT / 项目页 (BioProject) 说明。
- 追文献只对 "有论文标题或已发表、但 Assembly-Data 缺失" 的记录有价值 (老的 Sanger 之后、2012–2018 的 NGS 线粒体记录); 需在 annelid 等真实数据集上量化有多少条。

## 建议步骤
1. 在一个真实数据集(如 Annelida 线粒体 header)上统计: 无方法 且 有标题的记录数、去重后文献数。
2. 先做 "文献清单" (用户可自查, 部分文献难获取): 每篇一行 = 标题 / 作者 / 期刊 / 年 / DOI(Crossref, 复用 g2t-cite) / 受影响登录号列表 / 是否开放获取(Unpaywall) / PubMed 链接 / Google Scholar 检索链接(标题生成) → 输出 literature_to_check.xlsx。
3. 可选自动解析: 仅对开放获取全文 (PMC / Unpaywall PDF) 抽 Methods 段落, 用 mito_methods.yaml 规则匹配组装软件和平台, 结果写入新列 (assembler_basis = "paper full text (DOI)"), 全部标 "需人工确认", 不覆盖 GenBank 原值。
4. 不做: 自动爬 Google 搜索 (会被封, 违反条款); 付费墙文献不抓取。

## 待用户决定
- 步骤 2 先做、步骤 3 作为选项, 是否同意? 新增命令名建议 `g2t-cite --methods-list`, 或新子命令。
