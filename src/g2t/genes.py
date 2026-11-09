"""Gene types of g2t and their compartments — the one list the steps share.

Synonyms live in ``data/gene_dict.yaml``; the Entrez clauses for downloading in ``ncbi_genes.py``.
tests/test_gene_dictionary.py checks that the three agree. To add a gene type: add it to the
yaml, to one tuple below (its position is its matrix column) and to ``ncbi_genes.GENES``.
"""

from __future__ import annotations

# order = default column order of the matrix (g2t-organize --gene_order)
MITO_TYPES: tuple[str, ...] = (
    "mtgenome",
    "coi",
    "16s",
    "12s",
    "cob",
    "cox2",
    "cox3",
    "nd1",
    "nd2",
    "nd3",
    "nd4",
    "nd4l",
    "nd5",
    "nd6",
    "atp6",
    "atp8",
)
NUCLEAR_TYPES: tuple[str, ...] = ("18s", "28s", "its1-its2", "18-28s", "ef-1", "h3")
PLASTID_TYPES: tuple[str, ...] = (
    "plastome",
    "rbcl",
    "matk",
    "trnl-trnf",
    "psba-trnh",
    "rpl16",
    "rps16",
    "ndhf",
    "rpoc1",
    "ycf1",
    "rpl32-trnl",
    "trns-trng",
    "atpf-atph",
    "petb-petd",
    "psbe-petl",
    "rbcl-atpb",
    "trnk",
    "rps4",
    "atpb",
)
GENE_ORDER: tuple[str, ...] = MITO_TYPES + NUCLEAR_TYPES + PLASTID_TYPES

# whole-genome types; their records are copied into the columns of their compartment
GENOME_TYPES: frozenset[str] = frozenset({"mtgenome", "plastome"})

# gene names that people write into vouchers ('COI_ZMMU_WS399', 'KUN 12345 matK'); reconcile and
# curate strip them before comparing vouchers. Covers every single-locus type above plus common
# spellings (the test checks the coverage).
VOUCHER_GENE_TOKENS = (
    r"COI|CO1|COX1|COII|COX2|COIII|COX3|CYTB|COB|ND4L|NAD4L|ND\d|NAD\d|ATP\d|12S|16S|18S|28S|5\.8S|"
    r"H3|ITS1?|ITS2|EF-?1A?|EF1ALPHA|RRNL|RRNS|LSU|SSU|MTDNA|MITO|"
    r"RBCL|MATK|NDHF|RPOC1|YCF1|RPL16|RPS16|RPS4|ATPB|TRNK"
)
