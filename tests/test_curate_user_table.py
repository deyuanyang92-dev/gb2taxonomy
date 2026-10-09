"""Updating the NCBI matrix (Excel A) from the user's own Excel (B) with arbitrary column names -> Excel C."""

import pandas as pd

from g2t.curate import apply_user_table, detect_columns


def matrix():
    return pd.DataFrame([
        {"specimen_key": "Pc_WS399", "organism": "Priapulus caudatus", "voucher_standardized": "ZMMU:MSU:WS399",
         "voucher_as_submitted": "COI_ZMMU_MSU_WS399; WS399", "coi": "ON792938.1", "28s": "ON793000.1",
         "lat_lon": "", "geo_loc_name": "", "collection_date": "", "Ref1Title": "Direct Submission"},
        {"specimen_key": "Pc_WS3020", "organism": "Priapulus caudatus", "voucher_standardized": "ZMMU:MSU:WS3020",
         "voucher_as_submitted": "COI_ZMMU_MSU_WS3020", "coi": "ON792926.1", "28s": "",
         "lat_lon": "66.55 N 33.10 E", "geo_loc_name": "Russia", "collection_date": "", "Ref1Title": ""},
    ])


def user_table():
    return pd.DataFrame([
        {"标本号": "WS399", "COI 登录号": "ON792938", "28S accession": "ON793000",
         "生物种拉丁名": "Priapulus tuberculatospinosus", "纬度": "66.5512", "经度": "33.1023",
         "采集日期": "20190612", "论文题目": "Population genetics of Priapulus", "Depth (m)": "12"},
        {"标本号": "WS3020", "COI 登录号": "", "28S accession": "",
         "生物种拉丁名": "", "纬度": "-12.5", "经度": "-45.25",
         "采集日期": "2019-06-12", "论文题目": "", "Depth (m)": ""},
    ])


def test_detect_columns():
    m = detect_columns(user_table())
    assert m["标本号"] == "voucher"
    assert m["COI 登录号"] == "accession" and m["28S accession"] == "accession"
    assert m["生物种拉丁名"] == "organism"
    assert m["纬度"] == "latitude" and m["经度"] == "longitude"
    assert m["采集日期"] == "collection_date" and m["论文题目"] == "Ref1Title"
    assert m["Depth (m)"] == "user:Depth (m)"


def test_apply_user_table():
    out, log, problems, mapping = apply_user_table(matrix(), user_table(), key_column="specimen_key")
    a = out.set_index("specimen_key")
    assert a.loc["Pc_WS399", "organism"] == "Priapulus tuberculatospinosus"
    assert a.loc["Pc_WS399", "lat_lon"] == "66.5512 N 33.1023 E"
    assert a.loc["Pc_WS399", "collection_date"] == "12-Jun-2019"
    assert a.loc["Pc_WS399", "Ref1Title"] == "Population genetics of Priapulus"
    assert a.loc["Pc_WS399", "user:Depth (m)"] == "12"
    assert a.loc["Pc_WS3020", "lat_lon"] == "12.5 S 45.25 W"
    assert a.loc["Pc_WS3020", "collection_date"] == "12-Jun-2019"
    assert problems.empty
    assert set(mapping.columns) >= {"user_column", "used_as"}


def test_manual_mapping_overrides():
    t = user_table().rename(columns={"生物种拉丁名": "Taxon X"})
    out, log, _, mapping = apply_user_table(matrix(), t, key_column="specimen_key", column_map={"Taxon X": "organism"})
    assert out.loc[0, "organism"] == "Priapulus tuberculatospinosus"


def test_accessions_pointing_to_two_specimens_rejected():
    t = pd.DataFrame([{"accession": "ON792938", "accession 2": "ON792926", "species": "X y"}])
    out, log, problems, _ = apply_user_table(matrix(), t, key_column="specimen_key")
    assert log.empty and "point to different" in problems["problem"].iloc[0]


def test_cli_excel_a_plus_b_to_c(tmp_path):
    from g2t.curate import main
    a = tmp_path / "A.xlsx"
    with pd.ExcelWriter(a) as w:
        pd.DataFrame({"x": [1]}).to_excel(w, sheet_name="QC", index=False)
        pd.DataFrame([{"specimen_key": "Pc_WS399", "organism": "Priapulus caudatus",
                       "voucher_standardized": "ZMMU:MSU:WS399", "coi": "ON792938.1", "lat_lon": ""}]
                     ).to_excel(w, sheet_name="Matrix", index=False)
    b = tmp_path / "B.xlsx"
    user = pd.DataFrame([{"编号": "WS399", "纬度": "66.55", "经度": "33.10", "备注": "re-identified"}])
    user.to_excel(b, index=False)
    c = tmp_path / "C.xlsx"
    main(["-m", str(a), "-u", str(b), "-o", str(c), "--map", "编号=voucher"])
    sheets = pd.read_excel(c, sheet_name=None, dtype=str)
    assert list(sheets) == ["Matrix", "Changes", "Problems", "Column mapping", "Matrix (GenBank)"]
    m = sheets["Matrix"]
    assert m.loc[0, "lat_lon"] == "66.55 N 33.10 E" and m.loc[0, "user:备注"] == "re-identified"
    assert sheets["Matrix (GenBank)"]["lat_lon"].isna().all()


def test_values_starting_with_equals_are_written_as_text(tmp_path):
    import openpyxl

    from g2t.curate import write_curated_workbook, write_template
    df = pd.DataFrame({"specimen_key": ["a"], "Ref1Title": ["=HYPERLINK(\"http://x\")"]})
    write_curated_workbook(str(tmp_path / "c.xlsx"), df, pd.DataFrame(), pd.DataFrame(), df)
    write_template(df, str(tmp_path / "t.xlsx"))
    for f in ("c.xlsx", "t.xlsx"):
        for ws in openpyxl.load_workbook(tmp_path / f).worksheets:
            assert all(c.data_type != "f" for row in ws.iter_rows() for c in row)
