import json

import pandas as pd

from g2t import genomeasm
from g2t.litmethods import Cache


def test_batched_cached_lookup(tmp_path, monkeypatch):
    calls = []

    def fake_get(url, retries=4):
        calls.append(url)
        if "ena" in url:
            return b"accession\tsample_accession\nGCA_1\tS1\nGCA_2\tS2\n"
        a1 = {"assembly_name": "x", "assembly_level": "Scaffold", "assembly_method": "various"}
        a2 = {"assembly_method": "Hifiasm", "sequencing_tech": "PacBio"}
        reports = [{"accession": "GCA_1.1", "assembly_info": a1}, {"accession": "GCA_2.1", "assembly_info": a2}]
        return json.dumps({"reports": reports}).encode()

    monkeypatch.setattr(genomeasm, "_get", fake_get)
    t = pd.DataFrame({"accession": ["A.1", "B.1", "C.1"], "biosample": ["S1", "S2", "S3"], "organism": ["x"] * 3})
    cache = Cache(tmp_path / "c.json")
    r = genomeasm.build(t, cache)
    assert len(calls) == 2  # one ENA batch + one NCBI batch for all samples
    assert r.fillna("").set_index("accession").genome_method.to_dict() == {"A.1": "", "B.1": "Hifiasm", "C.1": ""}
    assert "no genome assembly" in r.set_index("accession").loc["C.1", "note"]
    genomeasm.build(t, Cache(tmp_path / "c.json"))
    assert len(calls) == 2  # second run is served from the cache
