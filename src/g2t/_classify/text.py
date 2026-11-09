"""Text helpers: phrase matching, DEFINITION cleaning, lengths, negative context."""

from __future__ import annotations

import logging
import re
from typing import Any

from g2t._classify.genedict import _phrase_pattern

logger = logging.getLogger("g2t.classify")  # one logger name for the whole step


def smart_match(phrase: str, text: str) -> bool:
    """Whole-word, case-insensitive phrase match."""
    return _phrase_pattern(phrase.lower()).search(text) is not None


def any_match(word_list: list[str], text: str) -> bool:
    return any(smart_match(w, text) for w in word_list)


def clean_length(length_str: str) -> int:
    """'550 bp' / '550' / 550.0 / '1,234 bp' -> int; anything else -> 0.
    (Stripping all non-digits turned 550.0 into 5500.)"""
    m = re.search(r"\d[\d,]*(?:\.\d+)?", str(length_str))
    if not m:
        return 0
    try:
        return int(float(m.group(0).replace(",", "")))
    except ValueError:
        return 0


# Words that look like gene keywords but are not marker genes.
NEGATIVE_CONTEXT: dict[str, re.Pattern] = {
    "16s": re.compile(
        r"16s\s+r?rna\s+(?:\((?:[^()]|\([^()]*\))*\)[\s-]*)?"
        r"(?:methyl|methylase|methyltransferase)",
        re.I,
    ),
    "12s": re.compile(r"12s\s+r?rna\s+(?:methyl|methylase|methyltransferase)", re.I),
    "h3": re.compile(
        r"(?:histone\s+)?h3\s*(?:k\d+|lysine|\S*(?:demethylase|methyltransferase|acetyltransferase|"
        r"deacetylase|kinase|chaperone))",
        re.I,
    ),
    "ef-1": re.compile(r"elongation\s+factor[\s-]*1[\s-]*(?:beta|gamma|delta|b|g|d)\b", re.I),
}
_SPECIMEN_TOKEN = re.compile(
    r"\b(?:isolate|voucher|strain|clone|specimen|haplotype|culture|individual|sample)"
    r"\s+(?:voucher\s+)?[^\s,;]+",
    re.I,
)


def clean_definition(definition: Any) -> str:
    """Lower-cased DEFINITION for keyword matching: specimen identifiers ('isolate CO2', 'voucher COI-12')
    and the English word 'its' are removed first, so they are not read as gene names."""
    text = str(definition or "")
    text = _SPECIMEN_TOKEN.sub(" ", text)
    text = re.sub(r"\bits\b", " ", text)  # lower-case possessive only; 'ITS' is kept
    return re.sub(r"\s+", " ", text).lower().strip()


def _negative_hit(gene: str, definition: str) -> bool:
    pat = NEGATIVE_CONTEXT.get(gene)
    return bool(pat and pat.search(definition))
