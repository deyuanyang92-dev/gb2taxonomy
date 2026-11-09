"""MatchConfig: the gene sets, dictionary and length limits of one run."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from g2t._classify.genedict import (
    _DEFAULT_MITO_GENES,
    _DEFAULT_NUCLEAR_GENES,
    _DEFAULT_PLASTID_GENES,
    default_gene_dict,
)

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


# =========================
# Configuration dataclass (holds the gene sets of a run)
# =========================
@dataclass
class MatchConfig:
    if_remove_matched: str = "no"

    length2_mtgenes: str = "150:50000"
    length2_ntgenes: str = "150:50000"
    length2_cpgenes: str = "150:300000"  # plastomes are ~120-220 kb; plastid records bypass the global filter
    global_length_range: str = "150:50000"

    # per-type length range; a type not listed has none ("none")
    gene_length_ranges: dict[str, str] = field(
        default_factory=lambda: {
            "mtgenome": "3000:500000",
            "plastome": "20000:300000",
            "18-28s": "5000:15000",
        }
    )

    add_assignment_reasons: bool = True
    record_conflict_reasons: bool = True

    # gene sets of this run
    mito_genes: set[str] = field(default_factory=lambda: _DEFAULT_MITO_GENES.copy())
    nuclear_genes: set[str] = field(default_factory=lambda: _DEFAULT_NUCLEAR_GENES.copy())
    plastid_genes: set[str] = field(default_factory=lambda: _DEFAULT_PLASTID_GENES.copy())

    # gene -> lower-case synonyms; each config owns its copy (see load_path_dict to add entries)
    gene_dict: dict[str, list[str]] = field(default_factory=default_gene_dict)

    def get_gene_length_range(self, gene: str) -> str:
        return self.gene_length_ranges.get(gene.lower(), "none")
