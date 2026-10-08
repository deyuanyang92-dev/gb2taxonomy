# g2t — GenBank to Taxonomy

[![version](https://img.shields.io/badge/version-v0.01-blue)](CHANGELOG.md) [![license](https://img.shields.io/badge/license-MIT-green)](LICENSE)

g2t downloads the GenBank records of a taxon from NCBI and turns them into a **specimen × gene matrix** for multi-locus phylogenetics and taxonomy: it reads the metadata of every record, assigns each record to a marker gene (COI, 16S, 18S, 28S, ITS, …), links the sequences of one specimen through its voucher, and writes one row per specimen with the accession of each gene.

中文使用说明见 [USAGE.md](USAGE.md)。Changes: [CHANGELOG.md](CHANGELOG.md). Known issues: [BUGS.md](BUGS.md).

## Installation

```bash
git clone https://github.com/deyuanyang92-dev/gb2taxonomy.git
cd gb2taxonomy
pip install -e ".[yaml]"
g2t --help
```

Requirements: Python ≥ 3.9, pandas ≥ 1.5, Biopython ≥ 1.80; pyyaml for custom gene dictionaries. The warning `Optional dependencies not installed: ete3` is harmless.

## Quick start

```bash
# 0. See what NCBI holds for the taxon, then download marker records
g2t-download -t Priapulidae -o gb --dry-run
g2t-download -t Priapulidae -o gb                 # -> gb/markers/batch_0001.gb ...

# 1–4. Metadata -> gene types -> vouchers -> reconciliation -> matrix
g2t -i gb/markers -o out --stream
# final matrix: out/organized_genes/organized_species_voucher.csv
```

Example (Priapulidae, NCBI txid37891, October 2026): 95,726 nuccore records, of which 69,704 are WGS contigs and 23,053 mRNA. The default selection keeps 548 marker records → 496 assigned to 7 gene types → 323 specimens, 117 with ≥ 2 genes.

## Usage

### Step 0 — download (`g2t-download`)

| Option | Selection |
|---|---|
| *(default)* | "markers": all records except WGS contigs, mRNA and RefSeq (predicted models and NC_/NR_ copies of INSDC records) |
| `--all` | everything, including WGS/mRNA/RefSeq (can be very large) |
| `--mito` | mitochondrial records |
| `--mitogenome` | complete mitochondrial genomes (10–30 kb) |
| `--gene COI,18S,28S` | given genes (COI COII COIII CYTB ND1 ND2 ND4 ND5 ATP6 12S 16S 18S 28S 5.8S ITS H3 EF1A RPB2; see `src/g2t/ncbi_genes.py`) |
| `--minlen N --maxlen N` | length range (bp) |
| `--query '…'` | any extra Entrez clause, e.g. `'Russia[Country]'` |
| `--include-wgs` / `--include-mrna` / `--include-refseq` | add back an excluded class |
| `--dry-run` | only report taxon composition, records per gene and the query |

- `-t` accepts a taxon name or an NCBI taxid. Each selection is saved in its own sub-directory (`markers`, `mito`, `mitogenome`, `gene-COI_18S`, …) together with `accessions.tsv` and `manifest.json` (query, date, counts).
- Records are fetched in batches of 200 by accession; every batch is checked for the expected number of records and retried on failure. Re-running the same command resumes.
- NCBI asks for an e-mail address: `-e you@example.org` or the variable `NCBI_EMAIL`. With an API key (`-k` or `NCBI_API_KEY`) NCBI allows 10 instead of 3 requests per second.
- `scripts/download_entrez.py` is a wrapper with the same options.

### Steps 1–4 — the pipeline (`g2t`)

```bash
g2t -i GB_DIR_OR_FILES -o OUT --stream [--resume] [--skip_reconcile] [--reconcile_min_confidence high]
```

| Step | Command | Input → output |
|---|---|---|
| 1 | `g2t-extract -i gb/markers -o s1` | GenBank → `final.csv` (one row per record: accession, length, definition, organism, voucher, isolate, locality, coordinates, date, references, …) |
| 2 | `g2t-classify -i s1/final.csv -o s2` | → `assigned_genes_types_all.csv` (gene type per record), `unmatched_sequences.csv`, `filtered_records.csv`, `record_status.csv` |
| 3 | `g2t-voucher -i s2/assigned_genes_types_all.csv -o s3` | → `updated_species_voucher.csv` (specimen key = organism + voucher) |
| 3b | `g2t-reconcile -i s3/updated_species_voucher.csv -o s3` | → `reconciled_species_voucher.csv`, `reconcile_report.csv` |
| 4 | `g2t-organize -i s3/reconciled_species_voucher.csv -o matrix.csv` | → specimen × gene matrix |

Output of the full pipeline:

```
OUT/
├── gb_metadata/final.csv                              step 1
├── labeled_genes/assigned_genes_types_all.csv         step 2
├── labeled_genes/unmatched_sequences.csv              step 2 (no gene type found)
├── labeled_genes/filtered_records.csv                 step 2 (records removed by a filter, with the reason)
├── labeled_genes/record_status.csv                    step 2 (every input record: assigned / unmatched / filtered + reason)
├── updated_species_vouchers/updated_species_voucher.csv     step 3
├── updated_species_vouchers/reconciled_species_voucher.csv  step 3b
├── updated_species_vouchers/reconcile_report.csv            step 3b (every candidate pair + decision)
└── organized_genes/organized_species_voucher.csv      step 4 (final matrix)
```

Gene types (step 2): `coi`, `cox2`, `cox3`, `cob`, `12s`, `16s`, `mtgenome`, `18s`, `28s`, `its1-its2`, `18-28s`, `ef-1`, `h3`. Synonyms are in `src/g2t/data/gene_dict.yaml`; add a gene or a spelling there, or pass your own file with `g2t-classify --path_dict`.

### Voucher reconciliation (step 3b)

Submitters often write the voucher of one specimen differently for each gene (`COI_ZMMU_MSU_WS399`, `28S_ZMMU_WS399`, `WS399`), so step 3 splits the specimen into one row per gene. Step 3b merges such groups only when the records support it:

- **Candidates**: same organism and same core identifier (gene names and collection acronyms removed; ≥ 4 characters). The same identifier in different organisms is reported, never merged.
- **Strong evidence**: the same publication (not "Direct Submission"); the same collection date; coordinates ≤ 0.01° apart.
- **Moderate evidence**: same first author, same collector, same detailed locality; coordinates ≤ 0.5° apart.
- **Conflict** (blocks the merge): incompatible dates, different countries, coordinates > 0.5° apart.
- **Decision**: ≥ 1 strong → merged (`high`); ≥ 2 moderate → merged (`medium`); otherwise not merged. If both groups already hold the same gene, only strong evidence merges them.

`--reconcile_min_confidence high` merges on strong evidence only; `--skip_reconcile` turns the step off. Merged rows carry `match_basis`, `match_confidence` and `match_evidence` into the matrix. In the Priapulidae data all 130 candidate pairs shared a publication and 94 specimens were merged; the identifier WS3020, used for *Halicryptus spinulosus* (28S) and *Priapulus caudatus* (COI, 16S), was reported and left unmerged.

### Python API

```python
import g2t
from g2t.download import DownloadOptions

r0 = g2t.download("Priapulidae", "gb", DownloadOptions(gene="COI,18S"))   # -> gb/gene-COI_18S/
res = g2t.run(input_files=[r0.out_dir], output_dir="out", stream=True)
print(res.success, res.final_output)

# single steps
g2t.extract(["gb/markers"], "s1")
g2t.classify("s1/final.csv", "s2")
g2t.voucher("s2/assigned_genes_types_all.csv", "s3")
g2t.reconcile("s3/updated_species_voucher.csv", "s3")
g2t.organize("s3/reconciled_species_voucher.csv", "matrix.csv")
```

## Limitations

1. **Gene type comes from the record's DEFINITION line, for the whole record.** A record that spans several regions gets one type. Multi-region rRNA records named "5.8S … ITS2 … 28S" or "18S … ITS1" are typed `its1-its2`, so their 18S or 28S part is missing from the 18s/28s columns (Priapulidae: 9 of 175 rRNA records, e.g. AY210840 with ~3.7 kb of 28S). See [BUGS.md](BUGS.md).
2. **Sequences are not cut into genes.** For a mitogenome (`mtgenome`) or an 18S–ITS–28S record (`18-28s`) the matrix copies the accession into every gene column it covers, but no per-gene sequence is extracted; this has to be done downstream before alignment.
3. **Length filter.** Records outside 150–50,000 bp (or without a length, or duplicated LocusIDs) are removed before classification. They are not lost silently: each is listed in `filtered_records.csv` with the reason, and `record_status.csv` gives every input record exactly one status. To classify them instead, widen the range with `g2t-classify … --length_range2_all 1:1000000`.
4. **Only 13 gene types.** Other loci (nuclear protein-coding genes, microsatellites, Hox genes, …) end up in `unmatched_sequences.csv` unless added to `gene_dict.yaml`.
5. **Vouchers.** Step 3 joins records of one organism whose voucher strings are identical, without further checks; two specimens that happen to share a code would be merged. Step 3b depends on the metadata submitted to GenBank: without a shared publication, date, coordinates or collector, variants of one voucher stay unmerged (it errs on the side of not merging).
6. **Names are taken as submitted.** Organism names are not checked against a taxonomic authority (WoRMS, NCBI synonyms); misidentified or outdated names stay as they are.
7. **Download.** NCBI nuccore only (no BOLD, ENA-only, or SRA data). `--gene` searches gene fields and title words, so it can miss records with unusual wording; a gene query also returns mitogenomes that contain the gene.
8. **Code status.** v0.01 is an early release. Tested on Python 3.13 with 232 unit tests; older modules still raise `ruff` style warnings.

## Citation

There is no paper on g2t yet. Please cite the software and version you used:

> Yang, D. (2026). *g2t: GenBank to Taxonomy* (version v0.01) [Computer software]. GitHub. https://github.com/deyuanyang92-dev/gb2taxonomy

```bibtex
@software{yang_g2t_2026,
  author  = {Yang, Deyuan},
  title   = {g2t: GenBank to Taxonomy},
  version = {v0.01},
  year    = {2026},
  url     = {https://github.com/deyuanyang92-dev/gb2taxonomy}
}
```

GitHub's "Cite this repository" button uses [CITATION.cff](CITATION.cff). Please also cite:

- **The sequence data**: the original publications of the GenBank records you use (columns `Ref1Authors`, `Ref1Title`, `Ref1Journal` of the output) and NCBI GenBank.
- **Biopython**, used to fetch and parse records: Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

## License

MIT — see [LICENSE](LICENSE).
