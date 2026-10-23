Chinese version: [USAGE.md](USAGE.md).

# g2t User Guide

## About

g2t (GenBank to Taxonomy) is a bioinformatics tool that arranges the molecular data of each specimen into a matrix for phylogenetic analysis.

In taxonomy and phylogenetics, several gene markers (e.g. COI, 16S, 18S) are often sequenced from the same specimen. g2t extracts this information from GenBank files, identifies the gene type of each record, links the sequences of one specimen through its voucher number, and produces a species × gene matrix.

## Requirements

**Core dependencies** (required):
- Python ≥ 3.9
- pandas ≥ 1.5
- biopython ≥ 1.80

**Optional dependencies**:
- pyyaml: reads the bundled `gene_dict.yaml`; without it an identical built-in copy is used
- openpyxl: reading and writing Excel files in Step 5 (`g2t-curate`)
- ete3: NCBI taxonomy lookup (for Class/Order/Family/Genus). If ete3 is not available (e.g. on Python 3.13), these columns are fetched once from NCBI Taxonomy by TaxonID and cached in `gb_metadata/taxonomy_cache.json`. Without network access the columns stay empty and a warning is printed.
- pytest: development tests

## Installation

### Option 1: use an existing environment (recommended)

If your environment already has pandas and biopython, install directly:

```bash
pip install -e /path/to/gb2taxonomy
```

### Option 2: create a new environment

```bash
# with mamba (recommended)
mamba create -n g2t python=3.10 -y
mamba activate g2t
pip install -e /path/to/gb2taxonomy

# or with conda
conda create -n g2t python=3.10 -y
conda activate g2t
pip install -e /path/to/gb2taxonomy

# or with pip + venv
python -m venv g2t_env
source g2t_env/bin/activate  # Linux/Mac
# g2t_env\Scripts\activate   # Windows
pip install -e /path/to/gb2taxonomy
```

### Option 3: install the dependencies only

```bash
pip install pandas biopython
```

### Check the installation

```bash
g2t --help
g2t --version
```

If a dependency is missing, the program says how to install it. Every command (`g2t`, `g2t-download`, `g2t-extract`, `g2t-classify`, `g2t-voucher`, `g2t-reconcile`, `g2t-organize`, `g2t-curate`) prints its version with `--version`.

---

## Quick start

```bash
# 0. See which records NCBI has for the taxon, then download the marker-gene records
g2t-download -t Priapulidae -o gb --dry-run
g2t-download -t Priapulidae -o gb                      # -> gb/markers/batch_0001.gb ...

# 1–4. Run the full pipeline (--stream is recommended for large files)
g2t -i gb/markers -o out --stream
# final matrix: out/organized_genes/organized_species_voucher.csv

# Use --resume to skip steps that are already done
g2t -i /path/to/files -o /path/to/output --resume
```

A small dataset to try offline is in [`examples/haplotaxidae`](examples/haplotaxidae): 347 GenBank records of Haplotaxidae (Annelida: Clitellata), about 4 s:

```bash
g2t -i examples/haplotaxidae/haplotaxidae.gb -o example_out
```

---

## Workflow

### Step 0: download

Downloads GenBank records of a taxon from NCBI. By default only "marker" records are downloaded: WGS contigs, mRNA and RefSeq records (predicted models, and NC_ and other copies that duplicate INSDC records) are excluded. In large taxa these classes often make up more than 95% of the records and are of no use for a multi-gene table.

```bash
g2t-download -t Priapulidae -o gb_out --dry-run                # taxon composition + records per gene + query
g2t-download -t Priapulidae -o gb_out                          # gb_out/markers/
g2t-download -t Priapulidae -o gb_out --mito                   # mitochondrial records
g2t-download -t Priapulidae -o gb_out --mitogenome             # complete mitochondrial genomes (≥10 kb, no upper limit; includes "genome assembly, organelle: mitochondrion")
g2t-download -t Priapulidae -o gb_out --gene COI,18S,28S       # selected genes (g2t/ncbi_genes.py)
g2t-download -t Priapulidae -o gb_out --gene COI --minlen 500
g2t-download -t Priapulidae -o gb_out --query 'Russia[Country]'
g2t-download -t Priapulidae -o gb_out --all                    # everything (can be very large)
```

- Records are downloaded in batches from a list of accessions (no WebEnv, which expires); failed requests are retried. Re-running the same command resumes the download: the accession list is fetched from NCBI again each time; a batch file is reused only if it holds exactly the accessions of that batch (new or withdrawn records are added or removed); batch files left over from a previous run are deleted.
- Each selection is stored in its own subdirectory, named after the options (markers / mito / mitogenome / gene-COI_18S / incl-wgs / all-mito / custom-<hash of the query> for `--query` …, or `--tag NAME`), together with `accessions.tsv` and `manifest.json` (query, date, checksum of the accession list, whether the download is complete). A directory made for another query is not reused (an error is raised; use another `--tag`). `--dry-run` writes no files.
- E-mail and API key: `-e/-k` or the environment variables `NCBI_EMAIL` / `NCBI_API_KEY` (optional; with a key, 10 requests per second).
- The old script `scripts/download_entrez.py` is now a wrapper around `g2t-download`. Its options are compatible, but the default selection (markers instead of all records), the output location (`<out>/<tag>/`) and the batch size (500) have changed; add `--all` for the old default selection.

### Step 1: extract metadata (extract)

Extracts record and specimen metadata from GenBank files (no sequences). Specimen fields are read from the `source` feature only; gene, rRNA and CDS features are not read.

```bash
g2t-extract -i /path/to/gb_files -o /path/to/out --stream
```

**Input**: GenBank flat files with the extension `.gb`, `.gbf`, `.gbk`, `.gbff` or `.genbank` (case-insensitive), or directories that contain them.

**Extracted fields**:
- Sequence: LocusID, ACCESSION, length, Definition
- Specimen: organism, specimen_voucher, isolate, strain, country, collected_by
- Taxonomy: TaxonID, Class, Order, Family, Genus
- Assembly: Assembly Method, Sequencing Technology

**Options**:
- `--stream`: streaming mode, for large files (recommended)
- `--batch`: process several files in parallel

**Output file**: `final.csv`

---

### Step 2: gene type (classify)

Assigns a gene type to each record from its Definition line.

```bash
g2t-classify -i final.csv -o /path/to/out
```

**Supported gene types**:

| Type | Gene markers |
|------|----------|
| Mitochondrial | coi, 16s, 12s, cob, mtgenome, cox3, cox2, nd1, nd2, nd3, nd4, nd4l, nd5, nd6, atp6, atp8 |
| Nuclear | 18-28s, its1-its2, ef-1, 28s, 18s, h3 |

A multi-region rRNA record whose Definition names both an ITS and 18S or 28S (e.g. "small subunit ribosomal RNA gene and internal transcribed spacer 1", "internal transcribed spacer 2 and large subunit ribosomal RNA gene") gets every type it covers, e.g. `gene_type` `18s,its1-its2` or `28s,its1-its2`; `Assignment_reason` notes "multi-region rDNA". A complete 18S–28S record stays `18-28s`.

**Output file**: `assigned_genes_types_all.csv`

---

### Step 3: specimen voucher (voucher)

Links the sequences of one specimen through its voucher number.

```bash
g2t-voucher -i assigned_genes_types_all.csv -o /path/to/out
```

**Voucher priority**:
1. specimen_voucher
2. isolate
3. culture_collection
4. clone
5. strain
6. ACCESSION (fallback)

**Output file**: `updated_species_voucher.csv`

---

### Step 3b: reconcile vouchers (reconcile)

Submitters often write the voucher of one specimen differently for different genes (`COI_ZMMU_MSU_WS399` / `28S_ZMMU_WS399` / `WS399`). Step 3 matches the strings literally and splits such a specimen into several rows. This step merges them based on evidence in the records:

| Evidence | Content |
|---|---|
| Candidates | A shared core identifier, taken from every identifier field (specimen_voucher, isolate, culture_collection, clone, strain), and either the same species or the same genus with an identical identifier (see below). Core: gene names, trailing notes such as "(holotype)" and punctuation are removed; year–number identifiers such as `2014-1234` are kept whole; in lists such as `WS399, WS400` each number counts separately; numbers with ≥5 digits also pair after their prefix is removed (`USNM 123456` ~ `123456`). Groups of the same genus whose identifiers are identical (ignoring case and punctuation, e.g. `SIO:BIC:A6302` = `SIO-BIC A6302`) are candidates even if the organism names differ (*Bathykurila* sp. / *Bathykurila* sp. AH-2023). Other pairs (different species, or a missing species name) are reported but not merged. An identical identifier in different genera is not merged and is flagged in `match_evidence` of both rows ("same voucher as … (different genus, not merged; check)"), with a warning |
| Strong | Identical identifier (ignoring case and punctuation; counts as moderate if both groups contain the same gene), the same publication (excluding Direct Submission), identical full collection date, coordinates ≤0.01° apart |
| Moderate | Same first author, same collector, same detailed locality, coordinates ≤0.5° apart, identical identifier when both groups contain the same gene |
| Conflict (blocks a merge) | The shared core is followed by different endings (lot or sub-sample numbers: `BIOUG06836-E01` / `-E02`, `SIO:BIC:A9392-1` / `-2`, `CE438` / `CE438_2`), dates differ, countries differ, coordinates >0.5° apart |

Decision: strong evidence → merged (high); ≥2 moderate items → merged (medium); otherwise not merged. If both groups already contain the same gene, only strong evidence that points to one specimen (same collection date or coordinates) can merge them; a shared publication is not enough, because one paper covers many specimens. Such pairs stay apart and are flagged in `match_evidence` ("similar voucher to … (both have the same gene, not merged; check)"). Before two specimens are merged, every pair of their records is checked. A chain A–B–C therefore does not join A and C if they conflict, differ in species (other than the genus rule above), or if they would put the same gene twice in one specimen without such evidence.

```bash
g2t-reconcile -i updated_species_voucher.csv -o /path/to/out [--min_confidence high]
```

Output: `reconciled_species_voucher.csv` (new columns species_voucher_g2t / voucher_core / match_basis / match_confidence / match_evidence) and `reconcile_report.csv` (evidence and decision for every candidate pair). The step is on by default in the full pipeline; `g2t --skip_reconcile` turns it off, and `--reconcile_min_confidence high` merges on strong evidence only.

---

### Step 4: build the matrix (organize)

Builds the species × gene matrix: one row per specimen, one column per gene marker.

```bash
g2t-organize -i updated_species_voucher.csv -o output.csv
```

**Output format**:

| species_voucher_new | organism | coi | 16s | 18s | ... |
|---------------------|----------|-----|-----|-----|-----|
| SpeciesA_voucher1 | Species A | LC123456 | LC123457 | LC123458 | |
| SpeciesB_voucher2 | Species B | LC123459 | | LC123460 | |

Each cell holds the LocusID of that specimen's sequence for that gene; several sequences are separated by semicolons. A record with several gene types (e.g. `18s,its1-its2`) appears in each of these columns.

**Output file**: `organized_species_voucher.csv`

---

### Vouchers in the matrix

After `organism` the matrix has three voucher columns:

- `voucher_as_submitted`: every specimen_voucher spelling of the specimen as it appears in GenBank, to show differences between submissions.
- `voucher_standardized`: one voucher per specimen, built from specimen_voucher only (empty for specimens that have only an isolate or other identifier). Gene-name prefixes are removed, the most complete spelling is used, and it is written in the INSDC `/specimen_voucher` format `[institution-code:[collection-code:]]specimen-id`. Examples: `COI_ZMMU_MSU_WS2585; 28S_ZMMU_WS2585; WS2585` → `ZMMU:MSU:WS2585`; `ZMMU:WS30980; ZMMU_WS30980` → `ZMMU:WS30980`.
- `voucher_note`: asks for a manual check when the spellings differ in more than prefix or separators.

### Step 5: curate metadata (curate, optional)

#### Purpose

Sequences are usually deposited in GenBank when a manuscript is submitted. Peer review often changes things afterwards: specimens are re-identified or given a new combination, species are described as new, the title changes, coordinates and dates are corrected. Only the submitter can update a GenBank record. INSDC welcomes "corrections of errors and update of the records by authors" ([INSDC policy](https://www.insdc.org/policy/)); NCBI asks submitters to send the publication data to GenBank once the paper is published ([GenBank overview](https://www.ncbi.nlm.nih.gov/genbank/)), and source information is updated by submitting a table ([Update GenBank records](https://www.ncbi.nlm.nih.gov/genbank/update/)). In practice many records are never updated, so the metadata in a downloaded dataset can be out of date or wrong. Misidentified and incompletely annotated public sequences are documented (Bridge et al. 2003, *New Phytologist* 160: 43–48, [doi:10.1046/j.1469-8137.2003.00861.x](https://doi.org/10.1046/j.1469-8137.2003.00861.x); Meiklejohn et al. 2019, *PLOS ONE* 14: e0217084, [doi:10.1371/journal.pone.0217084](https://doi.org/10.1371/journal.pone.0217084)).

Step 5 lets a taxonomist apply verified information to the specimen × gene matrix without changing GenBank and without losing the original values:

- **A**: the matrix built from GenBank, exactly as submitted (one row per specimen).
- **B**: your own table, in the format you already use. It holds the information you have verified: valid species names, published title and authors, corrected coordinates or dates, and extra fields such as water depth.
- **C**: the new matrix, A corrected with B, for downstream analysis.

Principles:

- **Per specimen**: each row of B corrects one specimen (one row of the matrix). The vouchers and accessions in B are used only to find that row.
- **Conservative**: if a row of B matches no specimen, matches several specimens, or its voucher and accession point to different specimens, it is not applied and the reason is listed.
- **Traceable**: C keeps the original values of A. Every change records the old value, the new value and the basis of the match; changed cells are highlighted in yellow. The change log is also the list of corrections to report to the submitter or to NCBI.
- **GenBank fields and formats**: values are written in GenBank format (`lat_lon` `66.55 N 33.10 E`, `collection_date` `12-Jun-2019`), so C can be compared directly with newly downloaded data. Darwin Core field names and export of an NCBI source-update table are planned.

g2t does not modify GenBank. C is your corrected copy; the records in INSDC remain the submitter's responsibility.

#### Usage

NCBI matrix (**A**) + your own table (**B**) → corrected matrix (**C**).

```bash
g2t-curate -m matrix_A.xlsx -u 我的表_B.xlsx -o 校正后_C.xlsx
g2t-curate -m matrix_A.xlsx -u B.xlsx -o C.xlsx --map "编号=voucher" --map "Lat=latitude"   # force a mapping when detection is wrong
```

(The file and column names in the first two commands are Chinese examples: "my table B", "corrected C", "number"; any language works.)

**B can be any table you already have** (Excel/csv, any column names, in Chinese or English; `--sheet` selects the worksheet). Columns are detected automatically:

| Column in B (examples) | Use |
|---|---|
| voucher, 凭证号, 标本号, catalog number | Find the specimen (any spelling) |
| Any accession column: `COI登录号`, `28S accession`, or a column whose values look like `ON792938.1` | Find the specimen; the accessions in one row must point to the same specimen |
| species, 拉丁名, 生物种拉丁名 (Latin name) | `organism` |
| Latitude + longitude (decimal) | Combined into GenBank `lat_lon` (`66.55 N 33.10 E`) |
| Collection date (`20190612`, `2019-06-12`) | `collection_date` in GenBank format (`12-Jun-2019`) |
| Locality, country, paper title, authors, journal, collector, identifier | `geo_loc_name`, `country`, `Ref1Title`, `Ref1Authors`, `Ref1Journal`, `collected_by`, `identified_by` |
| Any other column (e.g. `水深(m)`, water depth) | Added to C as a new column `user:<column name>` |

Worksheets of C: *Matrix* (changed cells in yellow), *Changes* (old value → new value, basis of the match), *Problems* (no match, several matching specimens, accession and voucher pointing to different specimens; none of these are applied), *Column mapping* (how each column of B is used), *Matrix (GenBank)* (A unchanged). A value equal to the one in A does not count as a change.

A pre-filled template can be used instead:

```bash
g2t-curate -m matrix.csv --template corrections.xlsx      # make a template: one row per specimen, current values filled in, plus a baseline worksheet
# edit the cells to be corrected in corrections.xlsx; empty cell = no change
g2t-curate -m matrix.csv -u corrections.xlsx -o curated.csv
```

- Locating a specimen: `accession` (an accession in any gene column; the version suffix can be omitted) and/or `voucher` (any spelling; e.g. `WS2585` finds `ZMMU_MSU_WS2585`). Add `organism_match` when the same voucher occurs in several species. All other columns are values to set (`organism`, `lat_lon`, `geo_loc_name`, `collection_date`, `Ref1Title`, `curation_source`, etc., or new columns).
- Output: `curated.csv` (with a `curated_fields` column), `curated_curation_log.csv` (every change: old value → new value, basis of the match), `curated_curation_problems.csv` (rows with no match, several matching specimens, or accession and voucher pointing to different specimens; none are applied). The original matrix is not changed.
- With a template, only cells you edited count as changes (compared with the `baseline (do not edit)` worksheet of the template). Cells you did not edit do not overwrite values that GenBank updated later. Cells that you edited and that GenBank has also changed are listed as conflicts and not applied.
- A cell cannot be emptied by curation (empty = no change).

---

## Output files

```
output/
├── gb_metadata/
│   └── final.csv                    # Step 1 output: metadata of all sequences
├── labeled_genes/
│   ├── assigned_genes_types_all.csv # Step 2 output: records with gene type
│   ├── unmatched_sequences.csv      # Step 2: records with no gene type
│   ├── filtered_records.csv         # Step 2: removed records and the reason (length, missing length, duplicate LocusID, etc.)
│   └── record_status.csv            # Step 2: one row per input record: assigned / unmatched / filtered + reason
├── updated_species_vouchers/
│   ├── updated_species_voucher.csv  # Step 3 output: records with specimen voucher
│   ├── reconciled_species_voucher.csv # Step 3b output: after voucher reconciliation
│   └── reconcile_report.csv         # Step 3b: evidence and decision for every candidate pair
├── organized_genes/
│   └── organized_species_voucher.csv # Step 4 output: final matrix
└── g2t_run.json                     # pipeline: date, command, versions (g2t, Python, Biopython, pandas, ete3, pyyaml), parameters, input files
```

---

## FAQ

### Q: How are sequences without a specimen_voucher handled?

The program tries isolate → culture_collection → clone → strain in turn, and finally uses the ACCESSION as a fallback.

### Q: How do I add a custom gene type?

Edit `src/g2t/data/gene_dict.yaml` (read only when pyyaml is installed; otherwise an identical built-in copy is used) and add the new gene type and its synonyms:

```yaml
gene_types:
  new_gene:
    synonyms:
      - synonym 1
      - synonym 2
```

Or use an external dictionary:

```bash
g2t-classify -i input.csv -o output --path_dict custom_dict.txt
```

`--path_dict` accepts `.txt` (one line per gene: `new_gene: synonym 1, synonym 2`), `.csv` (columns `gene,synonyms`, synonyms separated by commas) or `.py` (a variable `gene_dict`); `.yaml` is not accepted. The entries are added to the built-in dictionary.

A new gene type must also be listed in `g2t-classify --mtgenes_list` or `--ntgenes_list` (otherwise it is matched only in the second round, without a per-gene length range) and in `g2t-organize --gene_order` (otherwise it gets no column in the matrix). The `g2t` pipeline command does not pass these options, so run Steps 2 and 4 separately.

### Q: How do I run only some of the steps?

```bash
g2t -i input -o output --skip_extract    # skip Step 1
g2t -i input -o output --skip_classify   # skip Step 2
```

### Q: What if a dependency is reported missing?

The program detects missing dependencies and says how to install them. If you see:

```
ERROR: Missing required dependencies!
Missing: pandas, biopython
```

install them as suggested:

```bash
pip install pandas biopython
```

### Q: How do I record how a matrix was made?

The pipeline (`g2t`) writes `OUT/g2t_run.json` with the date, the command, the versions of g2t, Python, Biopython, pandas, ete3 and pyyaml, all parameters, and the input files (path and size). `g2t --version` prints the g2t version.

---

## Performance tips

- **Large files**: use `--stream` to avoid running out of memory.
- **Many files**: use `--batch` to process them in parallel.
- **Resuming**: use `--resume` to skip steps that are already done. If the reconcile settings differ from the previous run (stored in `OUT/.g2t_params.json`), Steps 3b and 4 are re-run. Re-running `g2t-classify` into the same directory first removes that step's previous outputs.

---

## Limitations and known issues

1. **Gene type is assigned from the whole DEFINITION line.** A multi-region rRNA record such as "18S … ITS1" or "ITS2 … 28S" is placed in every column it covers (`18s,its1-its2`, `28s,its1-its2`): the same accession appears in both columns and stands for the whole record, not for the separate region. Records that match other combinations of genes get all matched types (e.g. `16s,coi`, `28s,h3`) and `Conflict=True` (`Conflict_reason` "Multiple candidates…"); a mitochondrial record matching more than two genes and ≥3,000 bp long becomes `mtgenome`. See [BUGS.md](BUGS.md).
2. **Sequences are not cut into genes.** The accession of a mitochondrial genome (`mtgenome`) is copied into all mitochondrial gene columns (coi, 16s, 12s, cob, cox2, cox3), and that of an 18S–ITS–28S record (`18-28s`) into 18s, 28s and its1-its2, **whether or not the record actually annotates that gene** (e.g. FN689349 has no 12S annotation but appears in the 12s column). Restrict the columns with `g2t-organize --mtgenome_includes coi 16s`. No per-gene sequence is extracted; this has to be done separately before alignment.
3. **Length filter.** By default, records outside 150–50,000 bp (and records without a length or with a duplicate LocusID) are removed before classification. They are not lost silently: each is listed in `filtered_records.csv` with the reason, and `record_status.csv` gives every input record exactly one status. To classify them as well, widen all three ranges: `g2t-classify … --length_range2_all 1:1000000 --length2_mtgenes 1:1000000 --length2_ntgenes 1:1000000`.
4. **Only 22 gene types, matched by keywords.** Other loci (nuclear protein-coding genes, microsatellites, Hox genes, etc.) go to `unmatched_sequences.csv` unless added to the gene dictionary (see the FAQ on custom gene types). Specimen identifiers in the DEFINITION (`isolate CO2`, `voucher COI-12`) and the English word "its" are ignored, and several look-alikes are excluded (16S rRNA methyltransferase, histone H3 lysine/demethylase, elongation factor-1 beta/gamma). Keyword matching can still fail on unusual wording.
5. **Vouchers and metadata.** Step 3 merges records of one species whose voucher strings are identical, without further checks; two different specimens that happen to share a number are merged by mistake. Step 3b also merges identical identifiers (ignoring case and punctuation, in any identifier field) of one genus deposited under different names (*Bathykurila* sp. / *Bathykurila* sp. AH-2023); if both groups contain the same gene, only with the same date or coordinates; identical vouchers in different genera stay apart and are flagged in `match_evidence`. Step 3b depends on the metadata submitted to GenBank: without a shared publication, date, coordinates or collector, different spellings of one voucher are not merged (it errs on the side of not merging). When both groups already contain the same gene, a shared publication is not enough; the pair is flagged instead. Identifiers that differ only by a lot or sub-sample ending (`CE438` / `CE438_2`) are treated as different specimens. In the matrix, locality, date and other metadata are taken column by column from the first record of the specimen that has a value, so they may come from different records; `Conflict` is True if any record was flagged as a conflict.
6. **Species names are used as submitted.** They are not checked against WoRMS, NCBI synonyms or other taxonomic authorities; misidentifications and outdated names stay as they are unless corrected in Step 5.
7. **Download.** Only NCBI nuccore is searched (no BOLD, ENA-only data, or SRA). `--gene` searches gene fields and title keywords, so records with unusual wording can be missed; a gene query also returns mitochondrial genomes that contain the gene.
8. **Malformed files.** If Biopython cannot parse a record in a GenBank file, the records after it in that file are not read; the file is marked as partial in `extraction_report.json`, and a warning gives the expected and the parsed number of records.
9. **Curation (Step 5).** B only corrects specimens that are already in A. Specimens without a GenBank record are listed as "no match" and are not added to C. An empty cell means "no change", so B cannot empty a value in A. Column detection is a heuristic based on column names and on values that look like accessions; check the *Column mapping* sheet and use `--map` when detection is wrong. Values in B are taken as correct: the program does not check B against the literature or WoRMS.
10. **Code status.** v0.04 is an early release. 483 unit tests pass on Python 3.9, 3.12 and 3.13 (tested locally); `ruff` and `mypy` report no warnings.

---

## Citation

g2t has not been published in a paper yet. Please cite the software and the version you used:

> Yang, D. (2026). *g2t: GenBank to Taxonomy* (version v0.04) [Computer software]. GitHub. https://github.com/deyuanyang92-dev/gb2taxonomy

```bibtex
@software{yang_g2t_2026,
  author  = {Yang, Deyuan},
  title   = {g2t: GenBank to Taxonomy},
  version = {v0.04},
  year    = {2026},
  url     = {https://github.com/deyuanyang92-dev/gb2taxonomy}
}
```

"Cite this repository" on the right of the GitHub page reads [CITATION.cff](CITATION.cff). Please also cite:

- **Sequence data**: the original publications of the GenBank records you used (columns `Ref1Authors`, `Ref1Title`, `Ref1Journal` in the output) and NCBI GenBank.
- **Biopython** (used to download and parse the records): Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

---

## License

MIT License
