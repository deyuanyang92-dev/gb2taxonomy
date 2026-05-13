# NCBI GenBank Source Modifiers 格式参考

> 来源：https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html
> 更新页面：https://www.ncbi.nlm.nih.gov/genbank/update/
> 整理日期：2026-05-13

---

## 1. 两种使用场景的列名差异

| 场景 | 首列 | 地理位置 | 坐标 | 采集日期 |
|------|------|---------|------|---------|
| **BankIt 新提交** | `Sequence_ID` | `Country` | `Lat_Lon` | `Collection_date` |
| **记录更新 (gb-admin)** | `acc. num.` | `geo_loc_name` | `Lat_Lon` | `Collection_date` |

**注意**：更新已有记录时使用 `geo_loc_name`，不是 `Country` 或 `Country (geo_loc_name)`。

---

## 2. 更新已有记录的流程

1. 创建 Tab-delimited (.tsv) 文件，首列为 `acc. num.`
2. 仅包含需要更新的列（其他列留空或省略）
3. 发送邮件至 **gb-admin@ncbi.nlm.nih.gov**，附上 .tsv 文件
4. 必须是记录作者或获得授权

---

## 3. Source Modifiers 完整字段列表

### 必填字段（2024年12月起）

- **collection_date** — 采集日期
- **geo_loc_name** — 地理位置

### 常用字段及格式要求

| 字段 | 格式 | 示例 |
|------|------|------|
| `geo_loc_name` | `Country: Region, Locality` | `China: Sichuan, Chengdu` |
| `Lat_Lon` | `dd.dd N\|S ddd.dd E\|W` | `30.57 N 104.07 E` |
| `collection_date` | `DD-Mon-YYYY` / `Mon-YYYY` / `YYYY` | `08-Nov-2005`, `Dec-2000`, `1993` |
| `collected_by` | 自由文本 | `Smith, J.` |
| `specimen_voucher` | `institution-code:collection-code:specimen-id` | `SMNH:211581` |
| `isolation_source` | 自由文本 | `soil`, `freshwater` |
| `host` | 自由文本 | `Homo sapiens` |
| `isolate` | 自由文本 | `strain ABC123` |
| `altitude` | `value unit` | `1500 m` |
| `habitat` | 自由文本 | `marine` |
| `culture_collection` | `institution:accession` | `ATCC:12345` |
| `clone` | 自由文本 | `clone xyz1` |
| `strain` | 自由文本 | `strain K-12` |
| `note` | 自由文本 | 附加说明 |

### 其他字段

`identified_by`, `PCR_primers`, `sequenced_by`, `type_material`, `sex`, `breed`, `cultivar`,
`dev_stage`, `ecotype`, `focus`, `germline`, `haplogroup`, `haplotype`, `cell_line`,
`cell_type`, `tissue_type`, `mating_type`, `variety`, `sub_species`, `sub_strain`,
`serotype`, `biovar`, `biotype`, `group`, `pathovar`, `chemovar`, `sub_clone`,
`type_material`, `orig_filename`, `frequency`

---

## 4. 关键字段格式详解

### 4.1 geo_loc_name（地理位置）

格式：`Country: Region, Locality`

- 国家名使用 ISO 3166 标准英文名（首字母大写）
- 冒号 `:` 后跟省级区域，逗号分隔更细的地域
- 示例：
  - `Norway`（仅国家）
  - `Norway: Hordaland`（国家 + 省）
  - `Norway: Hordaland, Bergen, Mount Ulriken`（完整）

### 4.2 Lat_Lon（坐标）

格式：`latitude N|S longitude E|W`

- 十进制小数，保留 2-4 位
- 北纬用 `N`，南纬用 `S`；东经用 `E`，西经用 `W`
- 示例：`57.68 N 11.96 E`
- 不接受度分秒格式，需转换为十进制

### 4.3 collection_date（采集日期）

| 格式 | 适用场景 | 示例 |
|------|---------|------|
| `DD-Mon-YYYY` | 精确日期 | `08-Nov-2005` |
| `Mon-YYYY` | 仅知年月 | `Dec-2000` |
| `YYYY` | 仅知年份 | `1993` |
| `DD-Mon-YYYY to DD-Mon-YYYY` | 日期范围 | `01-Jan-2020 to 31-Dec-2020` |

月份缩写：Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec

### 4.4 specimen_voucher（标本凭证）

格式：`institution-code:collection-code:specimen-id`

- `specimen-id` 为必填部分
- `collection-code` 可选
- 示例：
  - `SMNH:211581`（机构:编号）
  - `USNM:Invertebrate Zoology:123456`（机构:部门:编号）
  - 个人收藏：`personal:John Smith:2020-001`

### 4.5 collected_by（采集人）

格式：`Lastname, Firstname` 或 `Lastname, F.`

- 多人用分号分隔：`Smith, J.; Doe, A.`
- 不缩写为 "et al."，应列出所有采集人

---

## 5. 已弃用字段（2025年1月起）

以下字段不再作为 Source Modifiers 使用：

Authority, Biotype, Biovar, Chemovar, Forma, Forma_specialis, Identified_by,
Pathovar, Pop_variant, Serogroup, Subclone, Subtype, Substrain, Type

---

## 6. 更新文件示例（Tab-delimited）

```
acc. num.	geo_loc_name	Lat_Lon	Collection_date	Collected_by	Specimen_voucher	Isolation_source
AB123456	Norway: Hordaland, Bergen	60.39 N 5.33 E	08-Nov-2005	Smith, J.	SMNH:211581	freshwater sediment
AB123457	China: Sichuan, Chengdu	30.57 N 104.07 E	1993	Wang, L.		personal:Li Wang:2020-001	soil
```

---

## 7. BankIt 新提交格式（参考对比）

```
Sequence_ID	organism	Country	Lat_Lon	Collection_date	Collected_by	Specimen_voucher	Isolation_source
seq1	Homo sapiens	USA: Maryland	38.95 N 77.08 W	01-Jan-2020	Smith, J.	personal:Smith:001	blood
```

**差异**：BankIt 首列为 `Sequence_ID`，列名用 `Country`；更新首列为 `acc. num.`，列名用 `geo_loc_name`。

---

## 8. 常见问题

**Q: 更新时是否需要包含所有字段？**
A: 不需要。只包含需要更新的字段即可，未列出的字段保持不变。

**Q: 空值如何处理？**
A: 留空即可（Tab 分隔但无内容）。NCBI 不会清除已有的值，除非明确要求。

**Q: 多条记录的更新能否放在同一文件？**
A: 可以。每行一条记录，首列为对应的 accession number。

**Q: 文件编码要求？**
A: 纯文本 Tab-delimited，推荐 UTF-8 编码。
