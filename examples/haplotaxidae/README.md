# Example: Haplotaxidae (347 GenBank records)

`haplotaxidae.gb` holds the 347 GenBank records (12S, 16S, 18S, 28S, COI) cited in Table S1 of
Martin P, Fend S, Martinsson S, Klinth M, Torii T, Erséus C (2024) Towards an integrative revision of
Haplotaxidae (Annelida: Clitellata). *Zoological Journal of the Linnean Society*.
[doi:10.1093/zoolinnean/zlae141](https://doi.org/10.1093/zoolinnean/zlae141).
The records were downloaded from NCBI nuccore by accession on 2026-05-13.

Run the pipeline (no network needed; without network or ete3 the Class/Order/Family/Genus columns stay empty):

```bash
g2t -i examples/haplotaxidae/haplotaxidae.gb -o example_out
```

Takes about 4 s. Expected: all 347 records assigned to a gene type; the matrix
`example_out/organized_genes/organized_species_voucher.csv` has one row per specimen, e.g. specimen CE1233
(*Haplotaxis* sp. SE) with 12S PP968574, 16S PP994580, 18S PP990634, 28S PP994523 and COI PP988405 in one row.
`tests/test_example_end_to_end.py` checks these results. The comparison of the whole matrix with Table S1 is
in `benchmark/`.
