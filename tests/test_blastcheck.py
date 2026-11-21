"""BLAST consensus: self hits dropped, independent units, identity gate, conflict needs margin + 2 units."""

import pandas as pd

from g2t import blastcheck as bc

LIN = {
    1: {"species": "Ophelia limacina", "genus": "Ophelia", "family": "Opheliidae", "order": "Opheliida",
        "class": "Polychaeta", "phylum": "Annelida"},
    2: {"species": "Ophelia bicornis", "genus": "Ophelia", "family": "Opheliidae", "order": "Opheliida",
        "class": "Polychaeta", "phylum": "Annelida"},
    3: {"species": "Notomastus latericeus", "genus": "Notomastus", "family": "Capitellidae", "order": "Opheliida",
        "class": "Polychaeta", "phylum": "Annelida"},
    4: {"species": "Capitella teleta", "genus": "Capitella", "family": "Capitellidae", "order": "Opheliida",
        "class": "Polychaeta", "phylum": "Annelida"},
}
DECL = {"species": "Ophelia neglecta", "genus": "Ophelia", "family": "Opheliidae", "order": "Opheliida",
        "class": "Polychaeta", "phylum": "Annelida"}


class FakeTax(bc.Taxonomy):
    def __init__(self):
        self.ok = True
        self.db = None

    def lineage(self, taxid):
        return LIN.get(taxid, {})


def hits(rows, q="HM746718.1"):
    return pd.DataFrame([{"qseqid": q, "sseqid": s, "pident": p, "length": 600, "qlen": 600, "slen": 600,
                          "qcovs": 100, "evalue": 0.0, "bitscore": p * 10, "staxids": str(t), "sscinames": ""}
                         for s, p, t in rows])


def run(rows, gene="coi", declared=DECL, q="HM746718.1"):
    return bc.check_query(q, hits(rows, q), declared, gene, FakeTax())


def test_self_hit_is_dropped():
    r = run([("HM746718.1", 100.0, 3), ("AB000001.1", 91.0, 1), ("AB000002.1", 90.5, 2)])
    assert r.best_pident == 100.0 or r.lca_taxon  # best_pident taken after dropping self
    assert r.best_pident < 100.0 and r.lca_taxon == "Ophelia"


def test_clean_record_has_no_conflict():
    r = run([("AB000001.1", 92.0, 1), ("AC000050.1", 91.8, 2)])
    assert r.gate == "actionable" and r.conflict_rank == "" and r.lca_taxon == "Ophelia"


def test_two_independent_foreign_units_is_a_family_conflict():
    r = run([("AB000001.1", 93.0, 3), ("AC009000.1", 92.8, 4), ("AD000001.1", 85.0, 1)])
    assert r.conflict_rank == "family" and r.lca_taxon == "Capitellidae"


def test_one_foreign_unit_alone_is_not_enough():
    r = run([("AB000001.1", 93.0, 3), ("AD000001.1", 85.0, 1)])
    assert r.conflict_rank == ""


def test_one_foreign_unit_from_own_batch_counts():
    r = run([("HM746701.1", 93.0, 3), ("AD000001.1", 85.0, 1)])
    assert r.same_batch_hit and r.conflict_rank == "family"


def test_batch_of_one_species_is_one_unit():
    rows = [(f"AB00{i:04d}.1", 93.0 - i * 0.01, 3) for i in range(1, 30)] + [("AD000001.1", 85.0, 1)]
    assert run(rows).conflict_rank == ""  # 29 sequences of one batch/species: one unit


def test_low_identity_is_no_close_reference_not_conflict():
    r = run([("AB000001.1", 70.0, 3), ("AC009000.1", 69.0, 4)])
    assert r.gate == "no_close_reference" and r.conflict_rank == ""


def test_declared_taxon_nearly_as_good_is_no_conflict():
    r = run([("AB000001.1", 91.5, 3), ("AC009000.1", 91.4, 4), ("AD000001.1", 90.5, 1)])
    assert r.conflict_rank == ""  # margin < 2


def test_conserved_marker_cannot_judge_genus():
    r = run([("AB000001.1", 98.0, 3), ("AC009000.1", 97.9, 4)], gene="18s")
    assert r.lca_rank in ("order", "class", "phylum")


def test_acc_key_and_batch():
    assert bc.acc_key("gb|HM746718.1|") == "HM746718"
    assert bc.same_batch("HM746718.1", "HM746701.1") and not bc.same_batch("HM746718.1", "AB746718.1")


def test_missing_rank_in_lineage_does_not_stop_the_consensus():
    """Annelid lineages have no 'order'; hits without a rank abstain instead of ending the walk."""
    saved = {k: dict(v) for k, v in LIN.items()}
    try:
        for v in LIN.values():
            v.pop("order")
        r = run([("AB000001.1", 93.0, 3), ("AC009000.1", 92.8, 4)])
        assert r.conflict_rank == "family" and r.lca_taxon == "Capitellidae"
    finally:
        LIN.update(saved)


def test_unclassified_hit_abstains():
    LIN[9] = {"phylum": "Annelida", "class": "Polychaeta"}
    try:
        r = run([("AB000001.1", 93.0, 3), ("AC009000.1", 92.8, 4), ("AE000001.1", 92.7, 9)])
        assert r.conflict_rank == "family"
    finally:
        del LIN[9]


def test_distant_hits_do_not_make_a_phylum_conflict():
    LIN[7] = {"phylum": "Mollusca", "class": "Gastropoda", "species": "Aplysia sp."}
    LIN[8] = {"phylum": "Mollusca", "class": "Bivalvia", "species": "Mytilus sp."}
    try:
        r = run([("AB000001.1", 78.0, 7), ("AC009000.1", 77.9, 8)])
        assert r.conflict_rank == "" and "distant_hit_only" in r.note
    finally:
        del LIN[7], LIN[8]


def test_read_hits_old_and_new_format(tmp_path):
    old = tmp_path / "old.tsv"
    old.write_text("q1\tgb|A1.1|\t99.0\t100\t100\t100\t100\t0.0\t200\t1\tfoo\n")
    new = tmp_path / "new.tsv"
    new.write_text("q1\tgb|A1.1|\t99.0\t100\t100\t100\t100\t0.0\t200\t1\tfoo\t5\t104\n")
    assert list(bc.read_hits(old)["qseqid"]) == ["q1"] and "sstart" in bc.read_hits(old)
    assert int(bc.read_hits(new)["send"][0]) == 104


def test_species_conflict_within_genus_needs_close_hits_of_another_species():
    # declared Ophelia neglecta; two independent units of O. limacina at 99%, none of the declared species
    r = run([("AB000001.1", 99.0, 1), ("AC009000.1", 98.9, 1), ("AD000001.1", 92.0, 2)])
    assert r.conflict_rank == "species" and r.lca_taxon == "Ophelia limacina"
    assert r.note.startswith("species_conflict")


def test_declared_species_as_good_as_the_foreign_one_is_no_species_conflict():
    lin6 = {**LIN[1], "species": "Ophelia neglecta"}
    tax = FakeTax()
    tax.lineage = lambda t: lin6 if t == 6 else LIN.get(t, {})
    r = bc.check_query("HM746718.1", hits([("AB000001.1", 99.0, 1), ("AC009000.1", 98.9, 1),
                                           ("AD000001.1", 99.0, 6)]), DECL, "coi", tax)
    assert r.conflict_rank == ""


def test_species_of_ignores_open_names():
    assert bc.species_of("Perinereis linea") == "Perinereis linea"
    for n in ("Perinereis sp.", "Nereis aff. zonata MALT-2022", "Nereididae environmental sample",
              "Nereis cf. zonata", "Nereididae"):
        assert bc.species_of(n) == ""


def test_discarded_hit_cannot_unlock_species_conflict():
    # two voting hits of O. limacina at 92%; a short 99% alignment of another species scores low and is outside TOP
    h = hits([("AB000001.1", 92.0, 1), ("AC009000.1", 92.0, 1), ("AD000001.1", 99.0, 2)])
    h.loc[2, "bitscore"] = 400.0
    r = bc.check_query("HM746718.1", h, DECL, "coi", FakeTax())
    assert r.best_pident == 92.0
    assert r.conflict_rank != "species"


def test_family_agreement_is_not_a_resolved_name():
    # declared species O. neglecta, only O. limacina references at 92%: genus agrees, species is unresolved
    r = run([("AB000001.1", 92.0, 1), ("AC009000.1", 91.8, 2)])
    assert r.agrees and r.confirmed_rank == "genus" and r.declared_rank == "species" and not r.resolved


def test_declared_genus_confirmed_is_resolved_and_conserved_markers_need_no_species():
    decl = {k: v for k, v in DECL.items() if k != "species"}
    assert run([("AB000001.1", 92.0, 1), ("AC009000.1", 91.8, 2)], declared=decl).resolved
    r = run([("AB000001.1", 99.5, 1), ("AC009000.1", 99.4, 2)], gene="18s")
    assert r.agrees and r.declared_rank == "family" and r.resolved


def _mini_taxdb(tmp_path):
    import sqlite3

    p = tmp_path / "taxa.sqlite"
    c = sqlite3.connect(p)
    c.execute("CREATE TABLE species (taxid INT, parent INT, spname TEXT, rank TEXT)")
    c.execute("CREATE TABLE synonym (taxid INT, spname TEXT)")
    c.execute("CREATE TABLE merged (taxid_old INT, taxid_new INT)")
    c.executemany("INSERT INTO species VALUES (?,?,?,?)", [
        (10, 1, "Nereididae", "family"), (11, 10, "Hediste", "genus"), (12, 11, "Hediste diversicolor", "species"),
        (13, 10, "Nereis", "genus"), (14, 13, "Nereis pelagica", "species"), (1, 1, "root", "no rank")])
    c.executemany("INSERT INTO synonym VALUES (?,?)", [(12, "Nereis diversicolor"), (12, "Nereis dual"),
                                                       (14, "Nereis dual")])
    c.execute("INSERT INTO merged VALUES (99, 12)")
    c.commit()
    return p


def test_accepted_name_resolves_synonym_and_merged_taxid(tmp_path):
    tax = bc.Taxonomy(_mini_taxdb(tmp_path))
    assert tax.accepted("Nereis diversicolor") == (12, "Hediste diversicolor")
    assert tax.accepted("Hediste diversicolor") == (12, "Hediste diversicolor")
    assert tax.accepted("Nereis dual") is None  # synonym of two taxa: ambiguous
    assert tax.accepted("Nope nope") is None
    assert tax.lineage(99).get("species") == "Hediste diversicolor"  # merged taxid


def test_family_only_units_do_not_give_a_species_consensus():
    # one species-labelled unit among nine family-only ones: genus / species stay unresolved (family kept)
    lin = dict(LIN)
    for t in range(10, 19):
        lin[t] = {"family": "Opheliidae", "order": "Opheliida", "class": "Polychaeta", "phylum": "Annelida"}

    class T(FakeTax):
        def lineage(self, taxid):
            return lin.get(taxid, {})

    rows = [("AB000001.1", 95.0, 1)] + [(f"AC0000{i:02d}.1", 94.9, 10 + i) for i in range(9)]
    r = bc.check_query("Q.1", hits(rows, "Q.1"), DECL, "coi", T())
    assert r.lca_rank == "family" and r.lca_taxon == "Opheliidae"


def test_single_foreign_hit_with_declared_taxon_absent_is_weak():
    r = run([("AB000001.1", 98.0, 3)], gene="coi")
    assert r.conflict_rank and r.weak and "declared" in r.note


def test_hits_lca_plain_band():
    import pandas as pd

    from g2t import blastcheck as bc

    class T:
        lin = {1: {"phylum": "Annelida", "family": "A", "genus": "G1", "species": "G1 x"},
               2: {"phylum": "Annelida", "family": "A", "genus": "G2", "species": "G2 y"},
               3: {"phylum": "Mollusca"}}

        def lineage(self, t):
            return self.lin[t]

    rows = [("q1", "AB1.1", 99, 1000, "1"), ("q1", "AB2.1", 98, 990, "2"), ("q1", "AB3.1", 80, 500, "3"),
            ("q1", "q1.1", 100, 1100, "3")]  # self hit dropped; 500 is outside the 2% band
    h = pd.DataFrame(rows, columns=["qseqid", "sseqid", "pident", "bitscore", "staxids"])
    r = bc.hits_lca(h, T()).iloc[0]
    assert (r.blast_lca_rank, r.blast_lca_taxon, r.blast_lca_n) == ("family", "A", 2)


def test_dataset_blast_flags_family_conflict(tmp_path):
    import random
    import shutil

    import pytest

    from g2t import blastcheck as bc

    if not shutil.which("blastn") or not shutil.which("makeblastdb"):
        pytest.skip("blast+ not installed")
    random.seed(1)
    base = "".join(random.choice("ACGT") for _ in range(400))

    def mut(s, n):
        s = list(s)
        for i in random.sample(range(len(s)), n):
            s[i] = random.choice([c for c in "ACGT" if c != s[i]])
        return "".join(s)

    seqs = {"A1.1": mut(base, 5), "A2.1": mut(base, 5), "A3.1": mut(base, 5)}
    fa = tmp_path / "coi.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    rows = {"A1": {"Family": "Foo", "Genus": "G", "Organism": "G a"},
            "A2": {"Family": "Foo", "Genus": "G", "Organism": "G b"},
            "A3": {"Family": "Bar", "Genus": "H", "Organism": "H c"}}
    d = bc.dataset_blast(fa, rows, threads=1).set_index("accession")
    assert bool(d.loc["A3.1", "ds_family_conflict"]) and d.loc["A3.1", "ds_family"] == "Foo"
    assert not d.loc["A1.1", "ds_family_conflict"]
