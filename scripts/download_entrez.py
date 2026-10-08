#!/usr/bin/env python3
"""
Download GenBank records by taxon (wrapper around ``g2t.download`` / ``g2t-download``).

Usage:
    # By taxon name (any rank) or NCBI taxid; marker records only (no WGS/mRNA/RefSeq)
    python scripts/download_entrez.py --taxon Hesionidae -e you@example.com -k API_KEY
    python scripts/download_entrez.py --taxon 6348 -o hesionidae_gb

    # Selections
    python scripts/download_entrez.py --taxon Priapulidae --mito
    python scripts/download_entrez.py --taxon Priapulidae --mitogenome
    python scripts/download_entrez.py --taxon Priapulidae --gene COI,18S,28S --minlen 300
    python scripts/download_entrez.py --taxon Priapulidae --dry-run      # counts only
    python scripts/download_entrez.py --taxon Priapulidae --all          # everything (can be huge)

Output:
    <output>/<tag>/batch_0001.gb, ... + accessions.tsv + manifest.json
    (tag = markers / mito / mitogenome / gene-COI_18S / ...). Re-running resumes.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from g2t.download import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
