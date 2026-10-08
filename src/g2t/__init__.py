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

__version__ = "1.1.0"


def run(*args, **kwargs):
    from g2t._pipeline import run as _run
    return _run(*args, **kwargs)


def download(*args, **kwargs):
    from g2t.download import download as _download
    return _download(*args, **kwargs)


def extract(*args, **kwargs):
    from g2t.extract import extract as _extract
    return _extract(*args, **kwargs)


def classify(*args, **kwargs):
    from g2t.classify import classify as _classify
    return _classify(*args, **kwargs)


def voucher(*args, **kwargs):
    from g2t.voucher import build_species_vouchers as _voucher
    return _voucher(*args, **kwargs)


def reconcile(*args, **kwargs):
    from g2t.reconcile import reconcile as _reconcile
    return _reconcile(*args, **kwargs)


def organize(*args, **kwargs):
    from g2t.organize import organize as _organize
    return _organize(*args, **kwargs)


__all__ = ["run", "download", "extract", "classify", "voucher", "reconcile", "organize"]
