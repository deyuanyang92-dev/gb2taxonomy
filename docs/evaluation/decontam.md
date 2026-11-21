# 污染 / 鉴定错误检查：BLAST 层、翻译检查、统一判定（2026-11-21）

## 设计
k-mer → 基因树 → BLAST 共识 → 翻译检查，合并为 `organized_genes/verdict.csv`（`g2t-screen --tree --blast --blast-db DB`）。
类别：likely_contaminant / likely_misidentified / suspect / numt_suspect / clean（需要 BLAST 正面证据且与声明谱系一致）/ unchecked。

## 真实数据（海蛹科 698 条，本地库 metazoa_marker，32 线程）
- 已知三条人工确认的错误全部被标出：OQ417162（k-mer + 树 + BLAST 三层，likely_misidentified）、HM746718、PX620130（BLAST 参考库里只有 1 个 80% 命中，仅基因树能看出；suspect）。
- LC670759–765（Travisia COI）六条翻译检查出内部终止子，标 numt_suspect，需人工确认是否为真 numt 或非 COI 片段。
- 各层互补：没有任何一层单独覆盖全部三条。

## 合成注入基准（T3，仅 BLAST + 翻译）
- 干净序列 200 条：BLAST 冲突 0 条。翻译层误报 3 条（1 条 COI，2 条 150 bp 片段）。
- numt（插入终止子 / 移码）检出 16/20。
- 发现并修复：单一近缘外科参照被漏（现 ≥97% 且比声明属高 ≥10 点即可报）；clean 未核对谱系（新增 agrees）；反向互补无标记（记录 sstrand）；ete3 库过旧致 no_hit（提示 ncbiupdate）。
- 修复后重测（同一数据）：外科混入 BLAST 冲突 17→29/30（suspect 29，clean 12→0，1 条因 taxid 缺失 unchecked）；标签互换 9→15/15；反向互补 0→20/20 标 `blast:reverse_complement`；numt 仍 16/20（4 条单碱基插入漏检）；干净序列 BLAST 冲突仍 0/200。代价：clean 需正面证据，COI 9/100、28S 21/50 变为 unchecked。同时修了 `read_hits` 读新格式文件时列名重复的 bug。
- Enchytraeidae 全量重跑：待做。

## 暂不做
rRNA HMM 本体检查（没有真值）；目标等级自动推断；study_group_id 批次合并（现为同前缀 500 号内）。
