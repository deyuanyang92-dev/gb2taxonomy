# scripts/

Helper scripts used while developing g2t. They are not installed with the package and are not needed to run it.
Scripts that read local data take its location from the environment variable `G2T_DATA_DIR` (default: the current directory).

| Script | Purpose | Status |
|---|---|---|
| `download_entrez.py` | download a taxon (thin wrapper around `g2t.download`) | use `g2t-download` |
| `download_requests.py`, `download_edirect.sh` | earlier download methods (requests / EDirect) | superseded by `g2t-download` |
| `download_by_accession.py` | download records from an accession list | kept: builds input for `benchmark/` from a published table |
| `compare_phyllodocida.py`, `validate_against_excel.py`, `compare_results.py` | first comparisons of g2t output with published tables | superseded by `benchmark/evaluate.py` |
| `merge_literature_metadata.py`, `generate_ncbi_update.py` | add metadata from a published table; write an NCBI source-update table | superseded by `g2t-curate` (Step 5) |
| `create_ncbi_template.py` | writes `docs/ncbi_submission_template.xlsx` | kept |
| `build_gene_dict_v2.py` | merged SynGenes synonym lists into a gene dictionary (needs a SynGenes checkout next to this repository) | kept for reference |
