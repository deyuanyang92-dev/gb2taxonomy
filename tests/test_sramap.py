import pandas as pd

from g2t import sramap
from g2t.litmethods import Cache


def test_build_record_runs_and_biosample_fallback(tmp_path, monkeypatch):
    calls = []

    def fake(q):
        calls.append(q)
        if "run_accession" in q:
            return [{"run": "SRR1", "strategy": "WGS", "platform": "ILLUMINA"}]
        return [{"run": "ERR9", "strategy": "Hi-C", "platform": "ILLUMINA"}]

    monkeypatch.setattr(sramap, "_query", fake)
    t = pd.DataFrame([
        {"accession": "A1", "sra_runs": "SRR1; SRR2", "biosample": ""},
        {"accession": "A2", "sra_runs": "", "biosample": "SAMEA1"},
        {"accession": "A3", "sra_runs": "", "biosample": ""},
    ])
    res = sramap.build(t, Cache(tmp_path / "c.json"))
    assert set(res.accession) == {"A1", "A2"}
    a1 = res[res.accession == "A1"].set_index("run")
    assert a1.loc["SRR1", "genomic"] == "yes" and a1.loc["SRR2", "note"] == "run not returned by ENA"
    a2 = res[res.accession == "A2"].iloc[0]
    assert a2.run_link == "biosample" and a2.genomic == "no"
