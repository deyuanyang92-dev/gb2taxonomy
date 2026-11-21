"""Translation check for protein-coding mitochondrial markers (numt / frameshift / wrong-gene screen).

A real mitochondrial COI / cob reads through in one frame with the invertebrate mitochondrial code
(NCBI table 5: TGA = Trp, ATA = Met, AGA / AGG = Ser; only TAA / TAG stop). Internal stops in every frame,
on both strands, mean a nuclear pseudogene (numt), a frameshifted / badly edited record, or not the gene.
The result is evidence for ``numt_suspect``, not a verdict.
"""

from __future__ import annotations

from dataclasses import dataclass

PROTEIN_CODING = {"coi", "cox1", "cob", "cytb", "coii", "coiii", "nd1", "nd2", "nd3", "nd4", "nd4l", "nd5", "nd6",
                  "atp6", "atp8"}
MIN_LEN = 150  # shorter fragments are not judged
_B = "TCAG"
_AA = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"  # table 5 (invertebrate mito)
CODON = {a + b + c: _AA[i] for i, (a, b, c) in enumerate((x, y, z) for x in _B for y in _B for z in _B)}
_COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


@dataclass
class TranslationResult:
    frame: int = 0  # 0..2 on ``strand``
    strand: str = "+"
    internal_stops: int = 0
    n_codons: int = 0
    flag: str = ""  # "" | numt_suspect | too_short | not_coding
    note: str = ""


def translate(seq: str, frame: int = 0) -> str:
    s = seq.upper().replace("U", "T")
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(frame, len(s) - 2, 3))


def check_translation(seq: str, gene: str) -> TranslationResult:
    """Best of six frames; internal stops (a terminal one is fine) decide."""
    res = TranslationResult()
    if gene.lower() not in PROTEIN_CODING:
        res.flag = "not_coding"
        return res
    seq = "".join(c for c in seq.upper() if c.isalpha())
    if len(seq) < MIN_LEN:
        res.flag = "too_short"
        return res
    best = None
    for strand, s in (("+", seq), ("-", seq.translate(_COMP)[::-1])):
        for f in range(3):
            p = translate(s, f).rstrip("*")
            stops = p.count("*")
            if best is None or stops < best[0]:
                best = (stops, f, strand, len(p))
    res.internal_stops, res.frame, res.strand, res.n_codons = best  # type: ignore[misc]
    if res.internal_stops:
        res.flag = "numt_suspect"
        res.note = f"{res.internal_stops} internal stop(s) in the best frame (table 5): numt, frameshift or wrong gene"
    return res
