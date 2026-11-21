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
