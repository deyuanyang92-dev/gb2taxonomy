#!/usr/bin/env python3
"""Convert published specimen tables into one long gold-standard table.

Output columns: dataset, row_id, species, voucher, specimen_level, gene, accession.

specimen_level is True when the row of the published table is one physical specimen
(it names a voucher/specimen number and cites no other study), so all its accessions
must end up in one matrix row. Rows that cite one or more publications may combine
sequences of different specimens; they are used for the gene-type check only.

Usage:
    python benchmark/build_gold.py --haplotaxis S1.xlsx --polynoidae TableS1.xlsx \
        --hesionidae supp.docx --opheliidae T0001.csv -o benchmark/gold/gold_standard.tsv
"""

from __future__ import annotations

import argparse
import re

import pandas as pd

ACC_RE = re.compile(r"^[A-Z]{1,2}_?\d{5,8}(\.\d+)?$")
GENE = {
    "12S": "12s",
    "16S": "16s",
    "16S rDNA": "16s",
    "18S": "18s",
    "18S rDNA": "18s",
    "28S": "28s",
    "28S rDNA": "28s",
    "28S sh/long": "28s",
    "COI": "coi",
    "COX1": "coi",
    "H3": "h3",
    "Histone 3": "h3",
    "CYTB": "cob",
    "Mitogenome": "mtgenome",
}


def _accessions(cell) -> list[str]:
    out = []
    for tok in re.split(r"[,;/\s]+", str(cell)):
        tok = tok.strip().rstrip("*")
        if ACC_RE.match(tok):
            out.append(tok.split(".")[0])
    return out


def _long(df: pd.DataFrame, dataset: str, species_col: str, voucher_col: str | None, level) -> list[dict]:
    rows = []
    for i, r in df.reset_index(drop=True).iterrows():
        voucher = str(r[voucher_col]).strip() if voucher_col else ""
        for col, gene in GENE.items():
            if col not in df.columns:
                continue
            for acc in _accessions(r[col]):
                rows.append(
                    {
                        "dataset": dataset,
                        "row_id": f"{dataset}:{i + 1}",
                        "species": str(r[species_col]).strip().rstrip(","),
                        "voucher": voucher,
                        "specimen_level": level(r),
                        "gene": gene,
                        "accession": acc,
                    }
                )
    return rows


def haplotaxis(path: str) -> list[dict]:
    df = pd.read_excel(path, header=3)
    df = df[df["Species"].notna() & df["Specimen no."].notna()]
    df["Voucher"] = df["Voucher"].astype(str) + " / " + df["Specimen no."].astype(str)
    return _long(df, "Haplotaxidae", "Species", "Voucher", lambda r: True)


def polynoidae(path: str) -> list[dict]:
    df = pd.read_excel(path, header=0)
    col = "Specimen Voucher or Source"

    def level(r):
        v = str(r[col])
        return bool(re.match(r"^[A-Z][A-Za-z-]+[- ][A-Z]*\d", v)) and ";" not in v

    return _long(df, "Polynoidae", "Species", col, level)


def hesionidae(path: str) -> list[dict]:
    import docx

    t = docx.Document(path).tables[0]
    data = [[c.text.strip() for c in row.cells] for row in t.rows]
    df = pd.DataFrame(data[1:], columns=data[0])
    # one cited study and no voucher: cannot be checked at specimen level
    return _long(df, "Hesionidae", "Species", "References", lambda r: False)


def opheliidae(path: str) -> list[dict]:
    df = pd.read_csv(path)
    return _long(df, "Opheliidae", "OTU", "Source/locality", lambda r: False)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--haplotaxis")
    p.add_argument("--polynoidae")
    p.add_argument("--hesionidae")
    p.add_argument("--opheliidae")
    p.add_argument("-o", "--output", required=True)
    a = p.parse_args()
    rows: list[dict] = []
    for fn, path in [
        (haplotaxis, a.haplotaxis),
        (polynoidae, a.polynoidae),
        (hesionidae, a.hesionidae),
        (opheliidae, a.opheliidae),
    ]:
        if path:
            rows += fn(path)
    out = pd.DataFrame(rows)
    # a complete mitogenome cited in several gene columns is one record: keep its mtgenome entry only
    mito = set(out.loc[out.gene == "mtgenome", ["row_id", "accession"]].itertuples(index=False, name=None))
    keep = [(g == "mtgenome") or ((r, a) not in mito) for r, a, g in zip(out.row_id, out.accession, out.gene)]
    out = out[keep]
    out.to_csv(a.output, sep="\t", index=False)
    print(out.groupby("dataset").agg(rows=("row_id", "nunique"), accessions=("accession", "nunique")))
    print("specimen-level rows:", out[out.specimen_level].row_id.nunique())


if __name__ == "__main__":
    main()
