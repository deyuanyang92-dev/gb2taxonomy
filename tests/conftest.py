import os

import pandas as pd
import pytest


@pytest.fixture
def sample_csv(tmp_path):
    """Create a minimal 3-row CSV with standard g2t columns."""
    csv_path = tmp_path / "sample.csv"
    df = pd.DataFrame(
        {
            "LocusID": ["LOC1", "LOC2", "LOC3"],
            "ACCESSION": ["ACC1", "ACC2", "ACC3"],
            "Definition": [
                "cytochrome oxidase subunit I",
                "16S ribosomal RNA",
                "18S ribosomal RNA",
            ],
            "Length": ["650", "550", "1800"],
            "MoleculeType": ["genomic DNA", "genomic DNA", "genomic DNA"],
            "organism": ["Species alpha", "Species beta", "Species gamma"],
            "gene_type": ["coi", "16s", "18s"],
        }
    )
    df.to_csv(csv_path, index=False)
    return str(csv_path)


@pytest.fixture
def sample_df(sample_csv):
    """Return the sample CSV as a DataFrame."""
    return pd.read_csv(sample_csv, dtype=str)
