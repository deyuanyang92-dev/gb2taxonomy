"""g2t-screen: k-mer distances and contamination / misidentification flags."""

import random

import pandas as pd

from g2t import screen as sc


def _mutate(s, p, rng):
    return "".join(c if rng.random() > p else rng.choice("ACGT".replace(c, "")) for c in s)


def test_kmers_are_strand_independent():
    a, b = sc.kmers("ACGTACGTTTGACCAGG", 7), sc.kmers("CCTGGTCAAACGTACGT", 7)
    assert (a == b).all() and sc.mash(a, b, 7) == 0.0
    assert sc.kmers("ACGTNNNNACGT", 5).size == 0  # windows with N are skipped


def test_mash_tracks_small_distances():
    rng = random.Random(1)
    s = "".join(rng.choice("ACGT") for _ in range(650))
    d = sc.mash(sc.kmers(s, 11), sc.kmers(_mutate(s, 0.02, rng), 11), 11)
    assert 0.01 < d < 0.03


def test_genus_of():
    assert sc.genus_of("Polyophthalmus pictus") == "Polyophthalmus"
    assert sc.genus_of("Opheliidae sp. A") == "" and sc.genus_of("uncultured annelid") == ""
    assert sc.genus_of("Armandia sp.", "Armandia") == "Armandia"


def _family(rng):
    """Two genera of 6 sequences each, ~15% apart; within a genus ~3%."""
    base = "".join(rng.choice("ACGT") for _ in range(650))
    a, b = _mutate(base, 0.08, rng), _mutate(base, 0.08, rng)
    ids, seqs, gen = [], [], []
    for g, root in (("Armandia", a), ("Polyophthalmus", b)):
        for i in range(6):
            ids.append(f"{g[:2].upper()}{i}.1")
            seqs.append(_mutate(root, 0.015, rng))
            gen.append(g)
    return ids, seqs, gen


def test_near_identical_to_another_genus_is_flagged():
    """Opheliidae COI: OQ417162 'Polyophthalmus pictus' is 0.4% from an Armandia sequence."""
    rng = random.Random(3)
    ids, seqs, gen = _family(rng)
    ids.append("OQ417162.1")
    seqs.append(_mutate(seqs[0], 0.004, rng))  # an Armandia sequence under a Polyophthalmus name
    gen.append("Polyophthalmus")
    rows = sc.screen_marker(ids, seqs, gen, fast=True)
    assert [(r["accession"], r["flag"], r["nearest_genus"]) for r in rows] == [
        ("OQ417162.1", "near_other_genus", "Armandia")]
    assert sc.screen_marker(ids, seqs, gen, fast=False) == []  # 18S / 28S: not judged by genus


def test_remote_sequence_is_flagged():
    rng = random.Random(4)
    ids, seqs, gen = _family(rng)
    ids.append("X1.1")
    seqs.append("".join(rng.choice("ACGT") for _ in range(650)))  # unrelated
    gen.append("Armandia")
    rows = sc.screen_marker(ids, seqs, gen, fast=True)
    assert [(r["accession"], r["flag"]) for r in rows] == [("X1.1", "remote")]


def test_matrix_gets_screen_flags(tmp_path):
    m = tmp_path / "m.csv"
    pd.DataFrame({"species_voucher_new": ["a", "b"], "coi": ["OQ417162.1", "AR0.1"], "16s": ["", ""]}).to_csv(
        m, index=False)
    flags = pd.DataFrame([{"gene": "coi", "accession": "OQ417162.1", "flag": "near_other_genus",
                           "note": "0.4% from Armandia"}])
    assert sc.mark_matrix(m, flags) == 1
    out = pd.read_csv(m, keep_default_na=False)
    assert out["screen_flags"].tolist() == ["OQ417162.1 coi: near_other_genus (0.4% from Armandia)", ""]


def test_fragment_against_complete_gene_is_not_penalised():
    rng = random.Random(5)
    full = "".join(rng.choice("ACGT") for _ in range(1231))
    frag = full[300:844]  # 544 bp of the same gene, identical
    a, b = sc.kmers(full, 11), sc.kmers(frag, 11)
    assert sc.pdist(a, b, 11) == 0.0 and sc.mash(a, b, 11) > 0.04  # Jaccard would call it 4.5% apart


def test_species_of():
    assert sc.species_of("Ophelia neglecta") == "Ophelia neglecta"
    for name in ("Ophelia sp. CMC01", "Ophelia aff. limacina CMC02", "Marionina argentea complex sp. CED",
                 "Enchytraeus environmental sample", "Opheliidae"):
        assert sc.species_of(name) == ""


def test_same_species_far_apart_is_flagged():
    """HM746718 'Ophelia neglecta' 16S is 24% from the 16S of the O. neglecta mitogenome OZ260669."""
    rng = random.Random(6)
    ids, seqs, gen = _family(rng)
    sp = [""] * len(ids)
    ids += ["OZ260669.1", "HM746718.1"]
    seqs += [_mutate(seqs[6], 0.03, rng), _mutate(seqs[6], 0.25, rng)]
    gen += ["Polyophthalmus", "Polyophthalmus"]
    sp += ["Polyophthalmus neglecta", "Polyophthalmus neglecta"]
    rows = sc.screen_marker(ids, seqs, gen, fast=True, species=sp)
    assert {r["accession"] for r in rows if r["flag"] == "same_species_far"} == {"OZ260669.1", "HM746718.1"}
    note = {r["accession"]: r["note"] for r in rows}
    assert "also named Polyophthalmus neglecta" in note["HM746718.1"]
