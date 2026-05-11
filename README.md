# g2t — GenBank to Taxonomy

A bioinformatics pipeline that processes GenBank sequence files into per-species summaries with gene type classifications and voucher identifiers.

## Installation

```bash
pip install .
# or for development
pip install -e .
```

**Dependencies:** Python >=3.9, pandas, biopython. Optional: pyyaml (for custom gene dictionaries).

## Quick Start

```bash
# Run full pipeline
g2t -i /path/to/genbank_files -o /path/to/output --stream

# Individual steps
g2t-extract -i /path/to/gb_files -o /path/to/out
g2t-classify -i final.csv -o /path/to/out
g2t-voucher -i assigned_genes_types_all.csv -o /path/to/out
g2t-organize -i updated_species_voucher.csv -o output.csv
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
| 1 | `g2t-extract` | Extract metadata from GenBank files (Biopython SeqIO) |
| 2 | `g2t-classify` | Classify gene types (13 types: COI, 16S, 18S, etc.) |
| 3 | `g2t-voucher` | Build species voucher identifiers |
| 4 | `g2t-organize` | Organize per-species gene summaries |

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
