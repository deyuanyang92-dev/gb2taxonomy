"""g2t-verify: rRNA regions from profile HMM hits, gene types the sequence supports, verdicts."""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from g2t import verify as vf
from g2t.verify import Hit


@pytest.mark.parametrize(
    "name, region",
    [("12S_rRNA", "12s"), ("16S_rRNA", "16s"), ("SSU_rRNA_eukarya", "18s"), ("18S_rRNA", "18s"),
     ("5_8S_rRNA", "5.8s"), ("LSU_rRNA_eukarya", "28s"), ("28S_rRNA", "28s"), ("5S_rRNA", None)],
)
def test_model_region(name, region):
    assert vf.model_region(name) == region


def test_its_from_5_8s_and_flanks():
    # 18S end 40, 5.8S 300-453, 28S from 700: ITS1 41-299, ITS2 454-699
    hits = [Hit("18s", 1, 40, 1792, 1831, 1831, 30), Hit("5.8s", 300, 453, 1, 154, 154, 90),
            Hit("28s", 700, 760, 1, 61, 3401, 35)]
    r = vf.regions(hits, 760)
    assert r.spans["its1"] == (41, 299, 0.0) and r.spans["its2"] == (454, 699, 0.0)
    assert vf.supported_types(r) == ["its1-its2"]  # 40 bp of 18S and 61 bp of 28S do not count


def test_overlapping_models_compete():
    """Mitochondrial 16S also aligns to the nuclear LSU model: the better model wins."""
    hits = [Hit("16s", 1, 480, 600, 1100, 1624, 400), Hit("28s", 50, 170, 900, 1020, 3401, 40)]
    r = vf.regions(hits, 480)
    assert set(r.spans) == {"16s"} and vf.supported_types(r) == ["16s"]


def test_full_cistron_is_18_28s():
    hits = [Hit("18s", 1, 1800, 30, 1831, 1831, 1500), Hit("5.8s", 2300, 2453, 1, 154, 154, 90),
            Hit("28s", 2800, 4000, 1, 1200, 3401, 900)]
    assert vf.supported_types(vf.regions(hits, 4000)) == ["18-28s"]


@pytest.mark.parametrize(
    "g2t, seq, any_hit, v",
    [
        (["12s"], ["16s"], True, "different"),  # Hesionidae MG523356-60
        (["28s", "its1-its2"], ["its1-its2"], True, "partly"),  # Enchytraeidae ITS records
        (["18s", "28s", "its1-its2"], ["18-28s"], True, "agree"),
        (["18-28s"], ["18-28s"], True, "agree"),
        (["28s", "its1-its2"], [], False, "unverifiable"),  # "ITS2 and 28S, partial", 237 bp
        (["16s"], [], False, "no_hit"),
        (["18s", "its1-its2"], [], True, "unverifiable"),  # KP745555: 48 bp of 18S, ITS1 without 5.8S
        ([], ["18s"], True, "new"),
        (["coi"], [], False, "agree"),
    ],
)
def test_verdict(g2t, seq, any_hit, v):
    assert vf.verdict(g2t, seq, any_hit) == v


def test_selection_takes_rrna_unmatched_and_conflicts():
    st = pd.DataFrame({"acc_base": ["A", "B", "C", "D", "E"],
                       "gene_type": ["coi", "16s", None, "28s,its1-its2", "coi"],
                       "status": ["assigned", "assigned", "unmatched", "assigned", "assigned"],
                       "Conflict": ["False", "False", "", "False", "True"]})
    assert vf.select(st) == {"B", "C", "D", "E"}
    assert vf.select(st, all_records=True) == {"A", "B", "C", "D", "E"}


# Stercutus niveus CE23080 (MF802217.1, Enchytraeidae): "ITS1, partial; 5.8S and ITS2, complete; LSU, partial".
# g2t types it 28s,its1-its2; the sequence holds ITS1-5.8S-ITS2 and only a short 28S end.
MF802217 = (
    "ACGTACGCACTCGGTCATCTCAGTGGGAATCGCTGGCCGCGCCGCCGGAGGTGAAACTCTGGTGCTGCGCGAGTCGGCGGGCGTCTGCTCTCGTAGGCGC"
    "TAAACAACTGGTGAACGCGTTCGGGTGCCTGGTCAGGCTCAAAGCATCGCACGGCTCTCGGCGATCGCCCGTCCGTAATTTTCTGTCTGCGACGGCGGAC"
    "GACCTGTATTCGCGGGTACCTATCTCGTCACGTCGCCGTGACGAGTGGAGCACGCTCTTGGCCAGGTTTAATGAATCGCCACCGGCGATCGCCCACTCGT"
    "CGTCGTCGTCAGACTTCTTTCACACCCCAACCACTACTCTGTTGTTACAAAACGCGAAAGCTCTCGGACGGACTGTCGCCGCCGTGTCGGCGCGTCCGTG"
    "TCGGGTCGATGCCCGTCCACAATTCGTTGCTCGGACGAGGCGGTGCGACCCGGCGACGCTGTCGAATGAAAATTGCAAAACGCAAATGATAACTCTAAGC"
    "GGTGGATCACTCGGCTCGTGTGTCGATGAAGAGCGCAGCCAGCTGCGTGAATTAATGTGAATTGCAGGACACATTGAACATCGATATCTTGAACGCATAT"
    "TGCGGCCTCGGGCATTCCCGAGGCCACGCCTGTCTCAGGGTCGGTTGAACGTCAATCGCGAGAAGCTCTCCGCTCGCGCATTGGGCGTCGCAGCATCGGG"
    "CGTCGGCGTCGGTCCTCTAGGATCGGCCGTCGAACGTCCGGCGGCTGCGTCGCCTCAAGCACAGACGGGTTGCTCGGTCGTCGGCGCGTCTGCCGAACCG"
    "GGACGGTGACTCGACTCTTGTAGAGCCCACCACCATCTCCCGGAGTTCGTGTCGGCTCGACTCGTCGGCGTACAGCGCAAGCGCCCACAAACGGCCGCTC"
    "GCGGCCAACAACACTGACTTCGACCTGAGATCAGACGAGATTACCCG"
)

REC = """LOCUS       MF802217                 947 bp    DNA     linear   INV 01-JAN-2018
DEFINITION  Stercutus niveus isolate CE23080 internal transcribed spacer 1, partial sequence; 5.8S
            ribosomal RNA gene and internal transcribed spacer 2, complete sequence; and large
            subunit ribosomal RNA gene, partial sequence.
ACCESSION   MF802217
VERSION     MF802217.1
SOURCE      Stercutus niveus
  ORGANISM  Stercutus niveus
            Eukaryota; Metazoa; Annelida; Clitellata; Enchytraeida; Enchytraeidae; Stercutus.
FEATURES             Location/Qualifiers
     source          1..947
                     /organism="Stercutus niveus"
ORIGIN
{origin}
//
"""


def _origin(seq):
    return "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10] for j in range(i, min(i + 60, len(seq)), 10))
                     for i in range(0, len(seq), 60))


@pytest.mark.skipif(importlib.util.find_spec("pyhmmer") is None, reason="pyhmmer not installed")
def test_end_to_end_with_bundled_rfam_models(tmp_path):
    gb = tmp_path / "gb"
    gb.mkdir()
    (gb / "batch_0001.gb").write_text(REC.format(origin=_origin(MF802217.lower())))
    out = tmp_path / "g2t" / "labeled_genes"
    out.mkdir(parents=True)
    pd.DataFrame({"ACCESSION": ["MF802217"], "status": ["assigned"], "gene_type": ["28s,its1-its2"]}).to_csv(
        out / "record_status.csv", index=False)
    t = vf.verify([str(gb)], str(tmp_path / "g2t"), str(tmp_path / "v"),
                  models=[vf.nuclear_models(tmp_path / "models")], cpus=1)
    r = t.iloc[0]
    assert r["sequence_type"] == "its1-its2" and r["verdict"] == "partly"
    assert "5.8s:" in r["regions"] and "its2:" in r["regions"]
    its2 = (tmp_path / "v" / "regions" / "its2.fa").read_text()
    assert its2.startswith(">MF802217.1 Stercutus niveus its2:")
    assert Path(tmp_path / "v" / "verify.csv").exists()


def test_short_record_of_one_region_counts():
    """GU453353 (175 bp, 28S D-region): the hit is short but it is the whole record."""
    r = vf.regions([Hit("28s", 3, 170, 400, 570, 3401, 27)], 175)
    assert vf.supported_types(r) == ["28s"]
    r = vf.regions([Hit("28s", 3, 70, 400, 470, 3401, 27)], 900)  # 68 bp of a 900 bp record: no
    assert vf.supported_types(r) == []


def test_stage1_models_follow_g2t_type():
    assert vf.stage1_regions(["16s"]) == {"5.8s", "16s"}
    assert vf.stage1_regions(["28s", "its1-its2"]) == {"5.8s", "18s", "28s"}
    assert vf.stage1_regions(["coi"]) is None and vf.stage1_regions([]) is None  # unknown: all models


def test_explained_needs_strong_hits_over_the_whole_record():
    assert vf.explained([Hit("28s", 1, 950, 1, 950, 3401, 300)], 960)  # 0.32 bits/bp, 10 bp left
    # LSU cross-hit on a mitochondrial 16S record (OP407585: 184 bp at 0.29 bits/bp of 552 bp): run the rest
    assert not vf.explained([Hit("28s", 200, 383, 900, 1080, 3401, 53)], 552)
    assert not vf.explained([Hit("28s", 1, 950, 1, 950, 3401, 250)], 960)  # weak: 0.26 bits/bp
    # ITS placed by a strong 5.8S counts as explained
    assert vf.explained([Hit("5.8s", 300, 453, 1, 154, 154, 90)], 700)


@pytest.mark.skipif(importlib.util.find_spec("pyhmmer") is None, reason="pyhmmer not installed")
def test_scan_dedups_and_caches(tmp_path, monkeypatch):
    models = vf.load_models([vf.nuclear_models(tmp_path / "models")])
    cache = vf.HitCache(tmp_path / "hits.sqlite")
    seqs = {"A.1": MF802217, "B.1": MF802217.lower()}  # identical sequences: scanned once
    calls = []
    real = vf._nhmmer
    monkeypatch.setattr(vf, "_nhmmer", lambda m, s, e: calls.append(len(s)) or real(m, s, e))
    first = vf.scan(seqs, models, cpus=2, cache=cache)
    assert sum(calls) == len(models) and first["A.1"] == first["B.1"]
    assert any(h.region == "5.8s" for h in first["A.1"])
    calls.clear()
    again = vf.scan({"C.1": MF802217}, models, cpus=2, cache=vf.HitCache(tmp_path / "hits.sqlite"))
    assert calls == [] and again["C.1"] == first["A.1"]  # from the cache, no nhmmer


def test_tile_hits_join_across_a_variable_gap():
    """PZ395572 (28S, 3,383 bp): tile hits 1-260, 211-419, then 738-893 after the D2 expansion segment."""
    hs = [Hit("28s", 1, 260, 44, 300, 3401, 110.9), Hit("28s", 211, 419, 254, 436, 3401, 147.1),
          Hit("28s", 738, 893, 636, 730, 3401, 54.7), Hit("28s", 902, 1206, 754, 1049, 3401, 254.1)]
    (m,) = vf.merge_tiles(hs)
    assert (m.start, m.end, m.hmm_from, m.hmm_to) == (1, 1206, 44, 1049)


def test_tiles_do_not_join_untiled_or_out_of_order_hits():
    two_58 = [Hit("5.8s", 432, 520, 1, 90, 154, 40), Hit("5.8s", 686, 839, 1, 154, 154, 90)]
    assert len(vf.merge_tiles(two_58)) == 2  # 5.8S is not tiled
    backwards = [Hit("28s", 1, 300, 2000, 2300, 3401, 200), Hit("28s", 400, 700, 100, 400, 3401, 200)]
    assert len(vf.merge_tiles(backwards)) == 2  # the model does not run on along the record


def test_minus_strand_tiles_join():
    hs = [Hit("16s", 1, 200, 1300, 1500, 1624, 120, "-"), Hit("16s", 180, 400, 1100, 1320, 1624, 130, "-")]
    (m,) = vf.merge_tiles(hs)
    assert (m.start, m.end, m.hmm_from, m.hmm_to) == (1, 400, 1100, 1500)


def test_ragged_ends_get_a_window_for_the_full_model():
    # KU728933 (1,255 bp): tiles found 16S 563-1206 (model 819-1532); the full model aligns from 55
    hs = [Hit("16s", 563, 1206, 819, 1532, 1624, 360)]
    assert vf.ragged(hs, 1255, {"16s"}) == [("16s", 1, 763), ("16s", 1006, 1255)]  # 49 bp left at the 3' end too
    # ends next to another region or at the model end are not ragged
    hs = [Hit("5.8s", 300, 453, 1, 154, 154, 90), Hit("28s", 454, 700, 1, 247, 3401, 120)]
    assert vf.ragged(hs, 700, {"28s"}) == []
    ext = vf.extend([Hit("16s", 563, 1206, 819, 1532, 1624, 360)], "16s", [Hit("16s", 55, 760, 400, 1010, 1624, 300)])
    assert (ext[0].start, ext[0].end, ext[0].hmm_from) == (55, 1206, 400)


def test_fast_reader_matches_biopython(tmp_path):
    from Bio import SeqIO

    f = tmp_path / "x.gb"
    f.write_text(REC.format(origin=_origin(MF802217.lower())) + REC.replace("MF802217", "MF802218").format(
        origin=_origin(MF802217[:120])))
    fast = list(vf.read_sequences([str(f)]))
    bio = [(r.id, r.annotations["organism"], str(r.seq).upper()) for r in SeqIO.parse(str(f), "genbank")]
    assert fast == bio and len(fast) == 2
    assert [a for a, _, _ in vf.read_sequences([str(f)], keep={"MF802218"})] == ["MF802218.1"]


def test_sequence_types_replace_annotation_types(tmp_path):
    d = tmp_path / "labeled_genes"
    d.mkdir()
    cols = {"ACCESSION": ["MG523356", "KP745545", "PX1", "AB1"], "Organism": ["a", "b", "c", "d"],
            "Assignment_reason": ["12s: keyword", "18s+its: keyword", "coi: keyword", ""],
            "gene_type": ["12s", "18s,its1-its2", "coi", ""]}
    a = pd.DataFrame(cols)
    a.iloc[:3].to_csv(d / "assigned_genes_types_all.csv", index=False)
    a.iloc[3:].to_csv(d / "unmatched_sequences.csv", index=False)
    pd.DataFrame({"ACCESSION": cols["ACCESSION"], "status": ["assigned"] * 3 + ["unmatched"],
                  "gene_type": cols["gene_type"], "reason": [""] * 4}).to_csv(d / "record_status.csv", index=False)
    t = pd.DataFrame({"accession": ["MG523356.1", "KP745545.1", "AB1.1"], "sequence_type": ["16s", "its1-its2", "28s"],
                      "verdict": ["different", "partly", "new"]})
    for _ in range(2):  # applying twice gives the same result
        st = vf.apply_to_assignments(t, d)
    assert st == {"changed": 2, "added": 1, "removed": 0}
    got = pd.read_csv(d / "assigned_genes_types_all.csv", dtype=str, keep_default_na=False)
    assert dict(zip(got["ACCESSION"], got["gene_type"])) == {"MG523356": "16s", "KP745545": "its1-its2", "PX1": "coi",
                                                             "AB1": "28s"}
    row = got.set_index("ACCESSION").loc["MG523356"]
    assert (row["annotation_type"], row["verify_verdict"]) == ("12s", "different")
    assert "12s -> 16s" in row["Assignment_reason"]
    assert (d / "assigned_genes_types_all.annotation.csv").exists()
    status = pd.read_csv(d / "record_status.csv", dtype=str, keep_default_na=False).set_index("ACCESSION")
    assert status.at["AB1", "status"] == "assigned" and status.at["MG523356", "gene_type"] == "16s"
