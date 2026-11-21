"""Phylogenetic screen: unrooted split test, record vs taxon decisions, distance veto, synonym genera."""

import importlib.util
import random

import pytest

from g2t import phylo as ph
from g2t.phylo import Tip

# Ophelia (study) with one contaminant HM746718 inside Notomastus (Capitellidae, background)
NWK = ("((O1:0.1,O2:0.1)1.00:0.1,(O3:0.1,O4:0.1)0.95:0.1,"
       "((HM746718:0.1,(N1:0.1,N2:0.1)0.99:0.1)0.99:0.2,(C1:0.1,C2:0.1)1.00:0.1)0.97:0.3);")


def _tips(names):
    info = {"O1": ("Ophelia limacina", "Ophelia", "Opheliidae", True),
            "O2": ("Ophelia borealis", "Ophelia", "Opheliidae", True),
            "O3": ("Ophelia neglecta", "Ophelia", "Opheliidae", True),
            "O4": ("Ophelia bicornis", "Ophelia", "Opheliidae", True),
            "HM746718": ("Ophelia neglecta", "Ophelia", "Opheliidae", True),
            "N1": ("Notomastus sp.", "Notomastus", "Capitellidae", False),
            "N2": ("Notomastus latericeus", "Notomastus", "Capitellidae", False),
            "C1": ("Capitella teleta", "Capitella", "Capitellidae", False),
            "C2": ("Capitella capitata", "Capitella", "Capitellidae", False)}
    return [Tip(n, *info[n]) for n in names]


def test_splits_and_unrooted_neighbours():
    names, parts = ph.splits(NWK)
    assert len(names) == 9 and all(0 <= s <= 1 for _, s in parts)
    i = names.index("HM746718")
    nb, sup = ph.neighbours(i, len(names), parts, (1 << len(names)) - 1)
    assert {names[j] for j in range(len(names)) if nb >> j & 1} == {"N1", "N2"} and sup == 0.99


def test_contaminant_flagged_outside_family():
    names, parts = ph.splits(NWK)
    rows, notes = ph.judge(_tips(names), parts, dict.fromkeys(names, 1.0), "Opheliidae")
    assert [(r["accession"], r["flag"]) for r in rows] == [("HM746718", "outside_family")]
    assert "4 other Ophelia sequences elsewhere" in rows[0]["note"] and notes == []


def test_low_coverage_and_low_support_are_not_judged():
    names, parts = ph.splits(NWK)
    cover = dict.fromkeys(names, 1.0)
    cover["HM746718"] = 0.2  # another region of the gene
    assert ph.judge(_tips(names), parts, cover, "Opheliidae")[0] == []
    assert ph.judge(_tips(names), parts, dict.fromkeys(names, 1.0), "Opheliidae", support=0.995)[0] == []


def test_whole_genus_elsewhere_is_a_taxon_note():
    """Travisia as a whole beside Capitellidae: a placement question for the genus, not record errors."""
    nwk = ("((O1:0.1,O2:0.1)1.00:0.1,((T1:0.1,T2:0.1)0.99:0.1,(C1:0.1,C2:0.1)0.98:0.1)0.97:0.2);")
    names, parts = ph.splits(nwk)
    info = {"O1": ("Ophelia limacina", "Ophelia", "Opheliidae", True),
            "O2": ("Ophelia borealis", "Ophelia", "Opheliidae", True),
            "T1": ("Travisia forbesii", "Travisia", "Opheliidae", True),
            "T2": ("Travisia pupa", "Travisia", "Opheliidae", True),
            "C1": ("Capitella teleta", "Capitella", "Capitellidae", False),
            "C2": ("Capitella capitata", "Capitella", "Capitellidae", False)}
    rows, notes = ph.judge([Tip(n, *info[n]) for n in names], parts, dict.fromkeys(names, 1.0), "Opheliidae")
    assert rows == []  # T1's neighbours include its congener T2


def test_genera_not_judged_on_conserved_markers():
    nwk = ("((O1:0.1,(P1:0.1,(A1:0.1,A2:0.1)0.99:0.1)0.99:0.1)0.95:0.1,"
           "(O2:0.1,(P2:0.1,(P3:0.1,P4:0.1)0.99:0.1)0.99:0.1)0.9:0.1);")
    names, parts = ph.splits(nwk)
    info = {n: ("Armandia sp." if n[0] == "A" else "Polyophthalmus pictus" if n[0] == "P" else "Ophelia sp.",
                {"A": "Armandia", "P": "Polyophthalmus", "O": "Ophelia"}[n[0]], "Opheliidae", True) for n in names}
    tips = [Tip(n, *info[n]) for n in names]
    flagged = [r["accession"] for r in ph.judge(tips, parts, dict.fromkeys(names, 1.0), "Opheliidae")[0]]
    assert flagged == ["P1"]  # OQ417162-like: Polyophthalmus inside Armandia, congener P2 elsewhere
    assert ph.judge(tips, parts, dict.fromkeys(names, 1.0), "Opheliidae", fast=False)[0] == []  # 18S / 28S


def _intermixed_rows(n_each: int) -> list[dict]:
    rows = []
    for i in range(n_each):
        rows.append({"gene": "16s", "accession": f"A{i}.1", "organism": "Ophiodromus sp.",
                     "flag": "inside_other_genus",
                     "note": "sits in Oxydromus (support 1.0); 5 other Ophiodromus sequences elsewhere",
                     "support": 1.0})
        rows.append({"gene": "16s", "accession": f"B{i}.1", "organism": "Oxydromus microantennatus",
                     "flag": "inside_other_genus", "note": "sits in Ophiodromus (support 0.91); 51 other ...",
                     "support": 0.91})
    rows.append({"gene": "16s", "accession": "DQ442575.1", "organism": "Leocrates chinensis",
                 "flag": "inside_other_genus", "note": "sits in Dalhousiella/Leocratides (support 0.93); ...",
                 "support": 0.93})
    return rows


def test_intermixed_genera_become_a_note():
    keep, notes = ph.intermixed_genera(_intermixed_rows(2), [])
    assert [r["accession"] for r in keep] == ["DQ442575.1"]
    assert {n["accession"] for n in notes} == {"A0.1", "A1.1", "B0.1", "B1.1"}


def test_one_record_each_way_may_be_a_swap_and_stays_flagged():
    keep, notes = ph.intermixed_genera(_intermixed_rows(1), [])
    assert {r["accession"] for r in keep} == {"A0.1", "B0.1", "DQ442575.1"}  # not hidden
    assert {n["accession"] for n in notes} == {"A0.1", "B0.1"}


def test_unresolved_tip_is_not_positive_evidence():
    names = ["P1", "P2", "O1", "O2", "O3"]
    tips = [Tip(n, "x", "G", "Opheliidae", True) for n in names[:2]] + \
           [Tip(n, "y", "H", "Other", False) for n in names[2:]]
    cover = dict.fromkeys(names, 1.0)
    assert ph.resolved_tips(tips, [], cover) == set()  # no supported split: unresolved


def _mut(s, p, rng):
    return "".join(c if rng.random() > p else rng.choice("ACGT".replace(c, "")) for c in s)


def test_close_congener_vetoes_a_tree_flag():
    """HQ025027 Travisia forbesii COI sat in Cirratulidae (1.0) but is 3.9% from another T. forbesii."""
    rng = random.Random(1)
    base = "".join(rng.choice("ACGT") for _ in range(600))
    trav, cirr = _mut(base, 0.04, rng), _mut(base, 0.35, rng)
    q = _mut(trav, 0.04, rng)
    entries = {"Q": (q, Tip("Q", "Travisia forbesii", "Travisia", "Opheliidae", True)),
               "T": (trav, Tip("T", "Travisia forbesii", "Travisia", "Opheliidae", True)),
               "C": (cirr, Tip("C", "Cirratulus sp.", "Cirratulus", "Cirratulidae", False))}
    rows = [{"accession": "Q", "others": ["C"], "note": "sits in Cirratulidae (support 1.0)"}]
    assert ph.distance_check(rows, {k: v for k, v in entries.items() if v[1].study}, entries) == []
    rows = [{"accession": "Q", "others": ["T"], "note": "sits in ..."}]  # closer to the group it sits in: kept
    assert len(ph.distance_check(rows, {"Q": entries["Q"]}, entries)) == 1


def test_orient_reverses_reversed_records():
    rng = random.Random(2)
    s = "".join(rng.choice("ACGT") for _ in range(500))
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    out, flipped = ph.orient({"a": s, "b": s[:400], "c": rc[50:450]})
    assert flipped == ["c"] and out["c"] in s


def test_orient_reversed_longest_reference_does_not_flip_the_batch():
    rng = random.Random(4)
    s = "".join(rng.choice("ACGT") for _ in range(900))
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    seqs = {"long_rc": rc, "b": s[:700], "c": s[100:800], "d": s[200:850], "e": s[50:600], "f": s[:400]}
    out, _ = ph.orient(seqs)
    assert all(v in s for v in out.values()) or all(v in rc for v in out.values())  # one common strand


@pytest.mark.skipif(not (ph.find_tool("famsa", "mafft") and ph.find_tool("FastTree", "fasttree", "iqtree3",
                                                                            "iqtree2", "iqtree")),
                    reason="FAMSA/MAFFT and FastTree/IQ-TREE not installed")
def test_build_tree_end_to_end(tmp_path):
    rng = random.Random(3)
    base = "".join(rng.choice("ACGT") for _ in range(500))
    seqs = {f"x{i}.1": _mut(base, 0.05 * (i % 3), rng) for i in range(8)}
    nwk, cover = ph.build_tree(seqs, ph.Tools.find(), tmp_path, threads=2)
    names, parts = ph.splits(nwk)
    assert sorted(names) == sorted(seqs) and min(cover.values()) > 0.9


def test_find_tool_missing():
    assert ph.find_tool("no-such-program-g2t") is None
    assert importlib.util.find_spec("Bio.Phylo") is not None


def test_resolved_needs_neighbourhood_that_agrees_with_name():
    # P1's supported neighbours are another family: settled, but it contradicts the name -> not confirmation
    names = ["P1", "P2", "O1", "O2", "O3"]
    tips = [Tip("P1", "x", "G", "Opheliidae", True), Tip("P2", "x", "H", "Opheliidae", True)] + \
           [Tip(n, "y", "K", "Other", False) for n in names[2:]]
    cover = dict.fromkeys(names, 1.0)
    parts = [(0b00111, 1.0)]  # P1, P2, O1 | O2, O3 ... P1's small side holds foreign O1
    assert "P1" not in ph.resolved_tips(tips, parts, cover)


def test_resolved_when_same_genus_neighbours():
    tips = [Tip(f"P{i}", "x", "G", "Opheliidae", True) for i in range(3)] + \
           [Tip(f"O{i}", "y", "K", "Other", False) for i in range(3)]
    cover = {t.name: 1.0 for t in tips}
    assert {"P0", "P1", "P2"} <= ph.resolved_tips(tips, [(0b000111, 1.0)], cover)


def test_representatives_keep_input_order():
    entries = {f"A{i:03d}": ("ACGT" * (i + 1), Tip(f"A{i:03d}", f"sp{i}", "G", "F", True)) for i in range(50)}
    out = ph.representatives(entries, cap=10, keep=set())
    assert list(out) == list(entries)


def test_support_scale_is_per_tree_not_per_node():
    # IQ-TREE UFBoot is 0-100: a node with 1 means 1 %, not 1.0
    names, parts = ph.splits("((a:1,b:1)1:1,(c:1,d:1)99:1,e:1);")
    assert sorted(s for _, s in parts) == [0.01, 0.99]
    names, parts = ph.splits("((a:1,b:1)0.3:1,(c:1,d:1)0.99:1,e:1);")
    assert sorted(s for _, s in parts) == [0.3, 0.99]


def test_resolve_support_defaults_by_tool():
    assert ph.resolve_support(None, "fasttree") == 0.90
    assert ph.resolve_support(None, "iqtree") == 0.95
    assert ph.resolve_support(95, "fasttree") == 0.95
    assert ph.resolve_support(0.8, "iqtree") == 0.8
