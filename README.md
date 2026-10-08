# g2t — GenBank to Taxonomy

A bioinformatics pipeline that downloads GenBank records of a taxon and turns them into a specimen × gene matrix, with gene type classification and evidence-based voucher reconciliation.

## Installation

```bash
pip install .
# or for development
pip install -e .
```

**Dependencies:** Python >=3.9, pandas, biopython. Optional: pyyaml (for custom gene dictionaries).

## Quick Start

```bash
# Step 0: download a taxon from NCBI (marker records by default; no WGS/mRNA/RefSeq)
g2t-download -t Priapulidae -o gb_out --dry-run           # composition + counts per gene
g2t-download -t Priapulidae -o gb_out                     # -> gb_out/markers/*.gb
g2t-download -t Priapulidae -o gb_out --mitogenome        # complete mitogenomes only
g2t-download -t Priapulidae -o gb_out --gene COI,18S --minlen 500

# Run full pipeline
g2t -i /path/to/genbank_files -o /path/to/output --stream

# Individual steps
g2t-extract -i /path/to/gb_files -o /path/to/out
g2t-classify -i final.csv -o /path/to/out
g2t-voucher -i assigned_genes_types_all.csv -o /path/to/out
g2t-reconcile -i updated_species_voucher.csv -o /path/to/out
g2t-organize -i reconciled_species_voucher.csv -o output.csv
```

## Python API

```python
import g2t

# Full pipeline
result = g2t.run(
    input_files=["/path/to/genbank_files"],
    output_dir="/path/to/output",
    stream=True,
)
print(f"Done: {result.steps_completed} steps, {result.elapsed:.1f}s")

# Individual steps
r1 = g2t.extract(["/path/to/gb_files"], "/path/to/out")
r2 = g2t.classify("/path/to/out/final.csv", "/path/to/out")
r3 = g2t.voucher("/path/to/out/assigned_genes_types_all.csv", "/path/to/out")
r4 = g2t.organize("/path/to/out/updated_species_voucher.csv", "output.csv")
```

## Pipeline Steps

| Step | CLI | Description |
|------|-----|-------------|
| 0 | `g2t-download` | Download a taxon from NCBI: markers / `--mito` / `--mitogenome` / `--gene` / length / custom query; resumable |
| 1 | `g2t-extract` | Extract metadata from GenBank files (Biopython SeqIO) |
| 2 | `g2t-classify` | Classify gene types (13 types: COI, 16S, 18S, etc.) |
| 3 | `g2t-voucher` | Build species voucher identifiers |
| 3b | `g2t-reconcile` | Merge voucher variants of one specimen (e.g. `COI_ZMMU_WS399` / `28S_ZMMU_WS399` / `WS399`) using evidence in the records |
| 4 | `g2t-organize` | Organize per-species gene summaries |

### Voucher reconciliation (step 3b)

Submitters often write the same voucher differently for each gene. Exact matching then splits one
specimen into one row per gene. `g2t-reconcile` merges such groups only when the records support it:

- **Candidates**: same organism and same core identifier (gene names / collection acronyms stripped, ≥ 4 characters). The same core in different organisms is reported, never merged.
- **Strong evidence**: shared publication title (not "Direct Submission"); identical collection date; coordinates within 0.01°.
- **Moderate evidence**: same first author, collector, or detailed locality; coordinates within 0.5°.
- **Conflict** (blocks the merge): incompatible dates, different countries, coordinates > 0.5° apart.
- **Decision**: ≥ 1 strong → merged (`high`); ≥ 2 moderate → merged (`medium`); otherwise not merged. If both groups already contain the same gene, only strong evidence can merge them. `--reconcile_min_confidence high` merges on strong evidence only; `--skip_reconcile` disables the step.

Every candidate pair and its decision is written to `updated_species_vouchers/reconcile_report.csv`; merged rows carry `match_basis`, `match_confidence` and `match_evidence` into the final matrix.

## Supported Gene Types

**Mitochondrial:** coi, 16s, 12s, cob, mtgenome, cox3, cox2
**Nuclear:** 18-28s, its1-its2, 28s, ef-1, 18s, h3

## CLI Reference

All commands support `-h` for full option lists. Key flags:

```bash
g2t -i INPUT -o OUTPUT [--stream] [--resume] [--quiet] [--skip_*] [--normalize_columns]
```

## License

MIT License. See [LICENSE](LICENSE).

## Citation

If you use g2t in your research, please cite:

```bibtex
@article{g2t2026,
  title = {g2t: a Python pipeline for GenBank-to-Taxonomy gene type classification and species voucher organization},
  author = {},
  journal = {Bioinformatics},
  year = {2026},
  doi = {}
}
```
