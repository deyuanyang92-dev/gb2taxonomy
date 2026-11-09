"""Package-level functions must stay callable after submodules are imported (review item A1)."""

import importlib

import pytest

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
    pd.DataFrame(
        [
            {
                "ACCESSION": "A",
                "LocusID": "A",
                "organism": "P c",
                "gene_type": "coi",
                "specimen_voucher": "WS399",
                "species_voucher_new": "k",
            }
        ]
    ).to_csv(f, index=False)
    assert g2t.reconcile(str(f), str(tmp_path)).success
    assert g2t.reconcile(str(f), str(tmp_path)).success


@pytest.mark.parametrize("name", ["g2t.extract", "g2t.classify"])
def test_step_modules_do_not_re_export_private_names(name):
    """Internals (underscore names) stay in g2t._extract / g2t._classify; the module keeps only its own helpers."""
    import ast
    import types
    from pathlib import Path

    mod = importlib.import_module(name)
    tree = ast.parse(Path(mod.__file__).read_text(encoding="utf-8"))
    own = {n.name for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    own |= {t.id for n in tree.body if isinstance(n, ast.Assign) for t in n.targets if isinstance(t, ast.Name)}
    leaked = [
        n
        for n, v in vars(mod).items()
        if n.startswith("_") and not n.startswith("__") and n not in own and not isinstance(v, types.ModuleType)
    ]
    assert not leaked, f"{name} re-exports private names: {leaked}"
    assert all(hasattr(mod, n) for n in mod.__all__)
