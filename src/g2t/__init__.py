"""
g2t — GenBank to Taxonomy
=========================
A pipeline: extract -> classify -> voucher -> reconcile -> organize.

Usage (CLI):
    g2t -i /path/to/gb_files -o /path/to/output --stream

Usage (Python API):
    import g2t
    g2t.run(input_files=["/path/to/gb"], output_dir="/path/to/out")
    result = g2t.extract(["/path/to/gb"], "/path/to/out")
"""

__version__ = "0.0.2"


# ``g2t.download``, ``g2t.extract`` ... are both submodules and the step functions. Loading a
# submodule rebinds ``g2t.<name>`` to the module object, so plain wrapper functions stop working
# after ``from g2t.download import DownloadOptions`` or ``g2t.run(...)``. Making the step modules
# callable keeps ``g2t.download(...)`` and ``import g2t.download as dl`` both valid.
import importlib as _importlib  # noqa: E402
import types as _types  # noqa: E402

_STEP_FUNCTIONS = {
    "download": "download", "extract": "extract", "classify": "classify",
    "voucher": "build_species_vouchers", "reconcile": "reconcile", "organize": "organize",
}


class _CallableModule(_types.ModuleType):
    def __call__(self, *args, **kwargs):
        return getattr(self, _STEP_FUNCTIONS[self.__name__.rsplit(".", 1)[1]])(*args, **kwargs)


for _name in _STEP_FUNCTIONS:
    _importlib.import_module(f"g2t.{_name}").__class__ = _CallableModule


def run(*args, **kwargs):
    from g2t._pipeline import run as _run
    return _run(*args, **kwargs)


__all__ = ["run", "download", "extract", "classify", "voucher", "reconcile", "organize"]
