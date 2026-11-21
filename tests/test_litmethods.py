from g2t import litmethods as lm

XML = ("<body><p>Illumina NovaSeq reads were assembled with GetOrganelle and NOVOPlasty.</p></body>"
       "<ref-list><ref>MitoZ: a toolkit (Meng 2019)</ref></ref-list>")


def test_methods_text_drops_reference_list(monkeypatch):
    monkeypatch.setattr(lm, "_get", lambda url, retries=3: XML.encode())
    txt = lm.methods_text("PMC1")
    assert "MitoZ" not in txt
    p = lm.parse_methods(txt)
    assert p["lit_seq_platforms"] == "Illumina"
    assert p["lit_assemblers"] == "GetOrganelle; NOVOPlasty"
