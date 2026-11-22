# 线粒体登录号的方法信息缺失：现状与提交建议

数据：EBP 线粒体集 6309 条（2026-11-22 运行，`g2t-methods`）。

## 现状

- 组装软件已知 4514 条（4506 条来自 GenBank 记录，8 条来自论文），**未知 1795 条（28%）**。
- 来源列 `assembler_source` / `platform_source` 逐值标明出处：`GenBank record` = 提交时已写；
  `paper PMID …` = 提交时没写、后来由论文补上；空 = 仍未知。
- 来自论文的只有 8 条，**这是检索覆盖不足，不是实际上限**：Europe PMC 仅覆盖约 5%；57 条论文提到登录号但未说明由其提交（按规则不采信）；
  23 条由 Methods 解析（需人工确认）；21 条全文获取失败。
- 1795 条未知记录多为批量提交，集中在少数单位：OceanOmics 531、Szeged 272、Trinity College Dublin 190、Surrey 149、Genoscope 72。
  提交年份 2018 年 394 条（可能已有论文）、2025–2026 年 965 条（论文大多还未发表）。
- 逐条清单：`unknown_assembler_1795.csv`（含提交日期、提交者、COMMENT、SRA 平台、核基因组方法）。

## 对后续提交的建议

1. 提交时就写 `##Assembly-Data-START##` 块：Assembly Method（软件名 + 版本）、Sequencing Technology。
2. COMMENT 里明确写线粒体基因组由哪个软件组装（如 MitoHiFi、OATK、GetOrganelle），不要只写核基因组流程。
3. 关联 BioSample 与 SRA run，保证 reads 可追溯。
4. 论文发表后，在 GenBank 记录里补 REFERENCE，论文里写明登录号与组装方法，便于后续检索。

## 说明

`g2t-methods` 只采信 GenBank 记录、SRA 和引用/报告了该登录号的论文，不做推断；缺失的值保持为空，不猜。
