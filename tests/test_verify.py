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
