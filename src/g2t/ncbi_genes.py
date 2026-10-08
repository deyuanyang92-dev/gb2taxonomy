"""Gene name -> Entrez query clauses, used by g2t.download (--gene).

Each gene searches both the [Gene] field and common title spellings (many older records
lack a /gene qualifier and name the gene only in DEFINITION).
To add a gene, add an entry to GENES; keys are case-insensitive.
"""

GENES = {
    # mitochondrial protein-coding
    "COI":   ['COI[Gene]', 'COX1[Gene]', 'CO1[Gene]', '"cytochrome c oxidase subunit I"[Title]',
              '"cytochrome oxidase subunit I"[Title]', '"cytochrome c oxidase subunit 1"[Title]',
              '"cytochrome oxidase subunit 1"[Title]'],
    "COII":  ['COII[Gene]', 'COX2[Gene]', '"cytochrome c oxidase subunit II"[Title]',
              '"cytochrome oxidase subunit 2"[Title]'],
    "COIII": ['COIII[Gene]', 'COX3[Gene]', '"cytochrome c oxidase subunit III"[Title]',
              '"cytochrome oxidase subunit 3"[Title]'],
    "CYTB":  ['CYTB[Gene]', 'COB[Gene]', '"cytochrome b"[Title]'],
    "ND1":   ['ND1[Gene]', 'NAD1[Gene]', '"NADH dehydrogenase subunit 1"[Title]'],
    "ND2":   ['ND2[Gene]', 'NAD2[Gene]', '"NADH dehydrogenase subunit 2"[Title]'],
    "ND4":   ['ND4[Gene]', 'NAD4[Gene]', '"NADH dehydrogenase subunit 4"[Title]'],
    "ND5":   ['ND5[Gene]', 'NAD5[Gene]', '"NADH dehydrogenase subunit 5"[Title]'],
    "ATP6":  ['ATP6[Gene]', '"ATP synthase F0 subunit 6"[Title]', '"ATP synthase subunit 6"[Title]'],
    # mitochondrial rRNA
    "12S":   ['12S[Title]', 'rrnS[Gene]', '"small subunit ribosomal RNA"[Title] AND mitochondrion[filter]'],
    "16S":   ['16S[Title]', 'rrnL[Gene]', '"large subunit ribosomal RNA"[Title] AND mitochondrion[filter]'],
    # nuclear
    "18S":   ['18S[Title]', '"18S ribosomal RNA"[Title]',
              '("small subunit ribosomal RNA"[Title] NOT mitochondrion[filter])'],
    "28S":   ['28S[Title]', '"28S ribosomal RNA"[Title]',
              '("large subunit ribosomal RNA"[Title] NOT mitochondrion[filter])'],
    "5.8S":  ['5.8S[Title]'],
    "ITS":   ['ITS1[Title]', 'ITS2[Title]', '"internal transcribed spacer"[Title]'],
    "H3":    ['H3[Gene]', '"histone H3"[Title]', '"histone 3"[Title]'],
    "EF1A":  ['EF1A[Gene]', 'EF1alpha[Gene]', '"elongation factor 1-alpha"[Title]',
              '"elongation factor 1 alpha"[Title]'],
    "RPB2":  ['RPB2[Gene]', '"RNA polymerase II"[Title]'],
}
ALIASES = {"COX1": "COI", "CO1": "COI", "COX2": "COII", "COX3": "COIII", "COB": "CYTB", "NAD1": "ND1",
           "RRNL": "16S", "RRNS": "12S", "SSU": "18S", "LSU": "28S", "EF1": "EF1A", "EF1ALPHA": "EF1A",
           "ITS1": "ITS", "ITS2": "ITS", "HISTONE3": "H3"}

MITO_GENES = ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND4", "ND5", "ATP6", "12S", "16S"]
NUCLEAR_GENES = ["18S", "28S", "5.8S", "ITS", "H3", "EF1A", "RPB2"]


def canon(name):
    k = name.strip().upper().replace("-", "").replace(" ", "")
    k = ALIASES.get(k, k)
    for g in GENES:
        if g.upper().replace("-", "").replace(".", "") == k.replace(".", ""):
            return g
    raise ValueError(f"Unknown gene '{name}'. Known: {', '.join(GENES)} (add it to g2t/ncbi_genes.py)")


def clause(gene):
    return "(" + " OR ".join(GENES[canon(gene)]) + ")"
