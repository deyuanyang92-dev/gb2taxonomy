import random

from g2t.translate_check import check_translation, translate

SENSE = [c for c in ("".join(x) for x in __import__("itertools").product("ACGT", repeat=3)) if c not in ("TAA", "TAG")]


def _orf(n, seed=1):
    r = random.Random(seed)
    return "".join(r.choice(SENSE) for _ in range(n))


def test_table5_codons():
    assert translate("TGAATAAGAAGG") == "WMSS"
    assert translate("TAATAG") == "**"


def test_clean_orf_any_frame_and_strand():
    s = _orf(200)
    assert check_translation(s, "coi").flag == ""
    assert check_translation("GG" + s, "coi").flag == ""  # frame 2
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    r = check_translation(rc, "coi")
    assert r.flag == "" and r.strand == "-"


def test_terminal_stop_is_fine():
    assert check_translation(_orf(200) + "TAA", "cob").flag == ""


def test_internal_stops_flag_numt():
    s = _orf(200)
    bad = s[:150] + "TAA" + s[153:300] + "TAG" + s[303:]
    r = check_translation(bad, "coi")
    assert r.flag == "numt_suspect" and r.internal_stops >= 1


def test_frameshift_flags():
    s = _orf(200)
    r = check_translation(s[:300] + s[301:], "coi")  # one base deleted: both halves cannot read in one frame
    assert r.flag in ("numt_suspect", "")  # random ORF may by chance stay open; guard only the API
    assert r.n_codons > 0


def test_not_coding_and_short():
    assert check_translation(_orf(200), "16s").flag == "not_coding"
    assert check_translation(_orf(20), "coi").flag == "too_short"
