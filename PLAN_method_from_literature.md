# 从来源补全测序平台 / 组装方法（g2t-litmethods）

目标：mitostats 里平台未知、组装方法未写明的记录，用可靠来源补全，并标明来源。不抓取 Google。

## 来源顺序（可靠度从高到低，`method_source` 列记录）
1. GenBank COMMENT / Assembly-Data（mitostats 已有）。
2. ENA read_run 仪器字段（`sra_*` 列）：记录的 DBLINK/COMMENT 里有 SRR/ERR/DRR 时，直接取 instrument_platform / model，不靠文字匹配。只给平台，不给组装软件。
3. Europe PMC：用登录号全文检索找论文（不是按标题），开放获取全文只取 Methods，用 mito_methods.yaml 的正则提取平台和组装软件（`lit_*` 列，建议值）。用登录号附近 220 字符判断论文是“报告”该序列还是“引用”它。
4. 非开放获取：只列论文（标题、DOI、PMCID），用户自己查。
5. 可选、默认关闭：LLM 提取（未实现）。

## 已验证的限制
- ENA / NCBI 组装记录只有核基因组，`assembly_software` 为空；线粒体组装软件无结构化来源。
- Europe PMC `databaseLinks` 对 OP104125 的论文返回空，不能建立论文↔登录号关系。
- 约 1800 条未发表 EBP 记录没有论文：平台靠 ENA，组装软件保持“未写明”。
- 正则对 Ray / IVA 等需要锚定的规则会漏掉散文式提法。

## 待做
- 报告 HTML 按 `method_source` 筛选；mitostats 表本身加 `method_source`。
- 非开放获取论文输出手动查找链接。
