import pandas as pd

from g2t.methodmerge import merge


def test_priority_record_then_sra_then_paper():
    st = pd.DataFrame([
        {"accession": "A", "organism": "x", "seq_platforms": "Illumina", "seq_category": "NGS", "assemblers": "GetOrganelle",
         "assembler_group": "GetOrganelle"},
        {"accession": "B", "organism": "y", "seq_platforms": "", "seq_category": "Unknown", "assemblers": "",
         "assembler_group": "No assembly info in record"},
        {"accession": "C", "organism": "z", "seq_platforms": "", "seq_category": "Unknown", "assemblers": "",
         "assembler_group": "No assembly info in record"},
    ])
    sra = pd.DataFrame([{"accession": "B", "genomic": "yes", "platform": "PACBIO_SMRT", "model": "Revio",
                         "run_link": "record"}])
    lit = pd.DataFrame([
        {"accession": "A", "lit_role": "reports", "lit_pmid": "1", "lit_assemblers": "NOVOPlasty", "lit_seq_platforms": ""},
        {"accession": "C", "lit_role": "reports", "lit_pmid": "2", "lit_assemblers": "MitoZ", "lit_seq_platforms": "Illumina"},
        {"accession": "B", "lit_role": "cited", "lit_pmid": "3", "lit_assemblers": "SPAdes", "lit_seq_platforms": ""},
    ])
    r = merge(st, sra, lit).set_index("accession")
    assert r.loc["A", "assembler"] == "GetOrganelle" and r.loc["A", "assembler_source"] == "GenBank record"
    assert r.loc["B", "platform"] == "PacBio" and r.loc["B", "assembler"] == ""
    assert r.loc["C", "assembler"] == "MitoZ" and "PMID 2" in r.loc["C", "assembler_source"]


def test_genome_columns_are_context_only():
    stats = pd.DataFrame({"accession": ["A1"], "organism": ["x"], "seq_platforms": [""], "seq_category": ["Unknown"],
                          "assemblers": [""], "assembler_group": [""]})
    gen = pd.DataFrame({"accession": ["A1", "A1"], "genome_assembly": ["GCA_1", "GCA_2"],
                        "genome_method": ["Hifiasm", "Hifiasm"], "genome_seq_tech": ["PacBio", "PacBio"]})
    r = merge(stats, genome=gen).iloc[0]
    assert r.genome_assembly == "GCA_1; GCA_2" and r.genome_method == "Hifiasm" and r.genome_seq_tech == "PacBio"
    assert r.assembler == ""
