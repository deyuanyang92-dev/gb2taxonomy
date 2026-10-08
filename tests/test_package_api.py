"""Package-level functions must stay callable after submodules are imported (review item A1)."""

import importlib

import g2t

NAMES = ["run", "download", "extract", "classify", "voucher", "reconcile", "organize"]


def test_functions_survive_submodule_imports():
    from g2t.download import DownloadOptions  # noqa: F401  (README example does this first)
    for n in NAMES:
        importlib.import_module(f"g2t.{n if n not in ('run',) else '_pipeline'}")
    for n in NAMES:
        assert callable(getattr(g2t, n)), n


def test_submodules_still_importable_as_modules():
    import g2t.download as dl
    assert hasattr(dl, "DownloadOptions") and callable(dl)


def test_reconcile_twice(tmp_path):
    import pandas as pd
    f = tmp_path / "v.csv"
    pd.DataFrame([{"ACCESSION": "A", "LocusID": "A", "organism": "P c", "gene_type": "coi",
                   "specimen_voucher": "WS399", "species_voucher_new": "k"}]).to_csv(f, index=False)
    assert g2t.reconcile(str(f), str(tmp_path)).success
    assert g2t.reconcile(str(f), str(tmp_path)).success
