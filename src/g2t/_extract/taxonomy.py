"""NCBI taxonomy: ete3 database, Entrez fallback and the lineage cache."""

from __future__ import annotations

import importlib.util
import json
import logging
import os
import re
import warnings

from Bio import Entrez

from g2t._extract.files import open_genbank

logger = logging.getLogger("g2t.extract")  # one logger name for the whole step


# ===============================
# Global state (encapsulated in TaxonomyService)
# ===============================
_LINEAGE_RANKS = ("class", "order", "family", "genus")
ENTREZ_TAXONOMY_BATCH = 200  # max TaxIDs per efetch request
TAXONOMY_CACHE_NAME = "taxonomy_cache.json"
_TAXON_XREF_RE = re.compile(r'/db_xref="taxon:(\d+)"')


class TaxonomyService:
    """Encapsulates NCBI taxonomy lookup state.

    Backends, in order of preference:
      1. ete3 (local NCBI taxonomy database) when it can be imported;
      2. otherwise Bio.Entrez ``efetch(db="taxonomy")`` in batches of at most
         200 TaxIDs.  Results are cached in ``<output_dir>/taxonomy_cache.json``.
    """

    def __init__(self):
        self.ete3_available: bool = False
        self.entrez_fallback: bool = False
        self.ncbi_taxa_cls = None
        self.cache: dict[str, dict[str, str]] = {}
        self.errors_reported: set = set()

    def update_database(self) -> bool:
        if not self.ete3_available:
            logger.warning("ete3 not installed")
            return False
        try:
            if self.ncbi_taxa_cls is None:
                from ete3 import NCBITaxa as C

                self.ncbi_taxa_cls = C
            ncbi = self.ncbi_taxa_cls()
            logger.info("Updating NCBI taxonomy database...")
            ncbi.update_taxonomy_database()
            logger.info("Update complete")
            return True
        except Exception as e:
            logger.error(f"Update error: {e}")
            return False

    def initialize(self, interactive: bool = True) -> bool:
        """Select a taxonomy backend.  Returns True when a backend is available.

        ``interactive`` is kept for backward compatibility; there is no prompt any
        more because the Entrez fallback makes "continue without taxonomy" moot.
        """
        self.ete3_available = False
        self.entrez_fallback = False
        reason = ""
        try:
            spec = importlib.util.find_spec("ete3")
        except (ImportError, ValueError):
            spec = None
        if spec is not None:
            try:
                from ete3 import NCBITaxa

                self.ncbi_taxa_cls = NCBITaxa
                self.ete3_available = True
                logger.info("ete3 module loaded")
                return True
            except Exception as e:  # ImportError / ModuleNotFoundError (e.g. no `cgi` on Python 3.13)
                reason = f"ete3 init error: {e}"
        else:
            reason = "ete3 not detected"
        self.entrez_fallback = True
        logger.warning(
            f"{reason}; falling back to NCBI Entrez (db=taxonomy, needs network) "
            f"for Class/Order/Family/Genus. Set NCBI_EMAIL / NCBI_API_KEY to identify yourself to NCBI."
        )
        return True

    def get_lineage(self, tax_id_str: str) -> dict[str, str]:
        if tax_id_str in self.cache:
            return self.cache[tax_id_str]
        if not self.ete3_available or not self.ncbi_taxa_cls:
            # Entrez fallback: lineages are bulk-loaded by prefetch_entrez(); never one request per record.
            return {}
        result, error = self._ete3_lineage(tax_id_str)
        if error and tax_id_str not in self.errors_reported:
            logger.warning(f"TaxonID '{tax_id_str}' query error: {error}")
            self.errors_reported.add(tax_id_str)
        self.cache[tax_id_str] = result
        return result

    def _ete3_lineage(self, tax_id_str: str) -> tuple[dict[str, str], str]:
        """Class/Order/Family/Genus from the local ete3 database, and the error text ("" when none)."""
        try:
            tax_id = int(tax_id_str)
        except ValueError:
            return {}, ""
        try:
            ncbi = self.ncbi_taxa_cls()
            lineage = ncbi.get_lineage(tax_id)
            if not lineage:
                return {}, ""
            # Not ncbi.get_rank()/get_taxid_translator(): ete3 3.1.3 quotes the ids with "...", which
            # SQLite >= 3.41 rejects ("no such column"), so every lineage came back empty.
            rows = ncbi.db.execute(
                f"SELECT taxid, rank, spname FROM species WHERE taxid IN ({','.join('?' * len(lineage))})",
                [int(t) for t in lineage],
            ).fetchall()
            info = {int(tid): (rank, name) for tid, rank, name in rows}
            result = {
                info[tid][0].capitalize(): info[tid][1]
                for tid in lineage
                if tid in info and info[tid][0] in _LINEAGE_RANKS
            }
            return result, ""
        except Exception as e:
            return {}, str(e)

    def prefetch_missing_local(self, tax_ids, cache_path: str | None = None) -> None:
        """ete3 backend: look up on NCBI Entrez the TaxonIDs the local database does not know (newer taxa)."""
        if not self.ete3_available or not self.ncbi_taxa_cls:
            return
        missing = []
        for t in tax_ids:
            t = str(t).strip()
            if t in self.cache:
                continue
            result, _ = self._ete3_lineage(t)
            if result:
                self.cache[t] = result
            else:
                missing.append(t)
        if missing:
            logger.info(f"Taxonomy: {len(missing)} TaxonID(s) not in the local ete3 database")
            self.prefetch_entrez(missing, cache_path, force=True)

    # ---- Entrez fallback -------------------------------------------
    @staticmethod
    def _load_disk_cache(cache_path: str | None) -> dict[str, dict[str, str]]:
        if not cache_path or not os.path.isfile(cache_path):
            return {}
        try:
            with open(cache_path, encoding="utf-8") as f:
                data = json.load(f)
            if not isinstance(data, dict):
                raise ValueError("not a JSON object")
            return {str(k): {str(a): str(b) for a, b in v.items()} for k, v in data.items() if isinstance(v, dict)}
        except Exception as e:
            logger.warning(f"Ignoring unreadable taxonomy cache {cache_path}: {e}")
            return {}

    @staticmethod
    def _save_disk_cache(cache_path: str, data: dict[str, dict[str, str]]) -> None:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
            tmp = cache_path + ".tmp"
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=1, sort_keys=True)
            os.replace(tmp, cache_path)
        except OSError as e:
            logger.warning(f"Could not write taxonomy cache {cache_path}: {e}")

    @staticmethod
    def _efetch_lineages(tax_ids: list[str]) -> dict[str, dict[str, str]]:
        """One Entrez request for ``tax_ids`` (at most ENTREZ_TAXONOMY_BATCH of them)."""
        email = os.environ.get("NCBI_EMAIL")
        api_key = os.environ.get("NCBI_API_KEY")
        if email:
            Entrez.email = email  # type: ignore[assignment]
        if api_key:
            Entrez.api_key = api_key  # type: ignore[assignment]
        saved = (Entrez.max_tries, Entrez.sleep_between_tries)
        Entrez.max_tries, Entrez.sleep_between_tries = 2, 2  # fail reasonably fast when offline
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)  # "Email address is not specified"
                handle = Entrez.efetch(db="taxonomy", id=",".join(tax_ids), retmode="xml")
                try:
                    records = Entrez.read(handle)
                finally:
                    handle.close()
        finally:
            Entrez.max_tries, Entrez.sleep_between_tries = saved
        out: dict[str, dict[str, str]] = {}
        for rec in records:
            nodes = list(rec.get("LineageEx", []))
            nodes.append({"Rank": rec.get("Rank", ""), "ScientificName": rec.get("ScientificName", "")})
            lineage = {}
            for node in nodes:  # ancestors first, the taxon itself last (same as ete3)
                rank = str(node.get("Rank", ""))
                if rank in _LINEAGE_RANKS:
                    lineage[rank.capitalize()] = str(node.get("ScientificName", ""))
            ids = [str(rec.get("TaxId", ""))] + [str(x) for x in rec.get("AkaTaxIds", [])]
            for tid in ids:
                if tid:
                    out[tid] = lineage
        return out

    def prefetch_entrez(self, tax_ids, cache_path: str | None = None, force: bool = False) -> None:
        """Bulk-load Class/Order/Family/Genus for ``tax_ids`` via Entrez (no-op unless fallback is active)."""
        if not self.entrez_fallback and not force:
            return
        wanted: list[str] = []
        seen = set()
        for t in tax_ids:
            t = str(t).strip()
            if t.isascii() and t.isdigit() and t not in seen:
                seen.add(t)
                wanted.append(t)
        disk = self._load_disk_cache(cache_path)
        for tid, lineage in disk.items():
            self.cache.setdefault(tid, lineage)
        missing = [t for t in wanted if t not in self.cache]
        if missing:
            logger.info(
                f"Taxonomy: requesting {len(missing)} TaxonID(s) from NCBI Entrez "
                f"({len(wanted) - len(missing)} already cached)"
            )
        for i in range(0, len(missing), ENTREZ_TAXONOMY_BATCH):
            chunk = missing[i : i + ENTREZ_TAXONOMY_BATCH]
            try:
                got = self._efetch_lineages(chunk)
            except Exception as e:
                pending = len(missing) - i
                logger.warning(
                    f"NCBI Entrez taxonomy lookup failed ({type(e).__name__}: {e}); "
                    f"Class/Order/Family/Genus will be empty for {pending} TaxonID(s). "
                    f"Re-run with network access to fill them."
                )
                break
            for tid in chunk:
                self.cache[tid] = got.get(tid, {})  # unknown ids -> {} (request succeeded)
        if cache_path:
            merged = dict(disk)
            merged.update({t: self.cache[t] for t in wanted if t in self.cache})
            if merged and (merged != disk or not os.path.isfile(cache_path)):
                self._save_disk_cache(cache_path, merged)


_taxonomy = TaxonomyService()


# ===============================
# NCBI taxonomy
# ===============================
def update_ncbi_taxonomy_database():
    return _taxonomy.update_database()


def initialize_ncbi_taxa(interactive: bool = True) -> bool:
    return _taxonomy.initialize(interactive)


def get_taxonomy_lineage(tax_id_str: str) -> dict[str, str]:
    return _taxonomy.get_lineage(tax_id_str)


def collect_taxon_ids(gb_files: list[str]) -> list[str]:
    """Unique ``/db_xref="taxon:<id>"`` values in the files (cheap line scan, no full parse)."""
    ids: dict[str, None] = {}
    for fp in gb_files:
        try:
            with open_genbank(fp) as f:
                for line in f:
                    if "taxon:" in line:
                        for tid in _TAXON_XREF_RE.findall(line):
                            ids[tid] = None
        except OSError:
            continue
    return list(ids)


def prefetch_taxonomy(gb_files: list[str], output_dir: str, include_taxonomy: bool = True) -> None:
    """Fill the taxonomy cache in bulk (Entrez fallback, or the TaxonIDs missing from the local ete3 database)."""
    if not include_taxonomy or not (_taxonomy.entrez_fallback or _taxonomy.ete3_available):
        return
    ids = collect_taxon_ids(gb_files)
    if not ids:
        return
    cache_path = os.path.join(output_dir, TAXONOMY_CACHE_NAME)
    if _taxonomy.entrez_fallback:
        _taxonomy.prefetch_entrez(ids, cache_path)
    else:
        _taxonomy.prefetch_missing_local(ids, cache_path)


def _init_worker(cache_dict, ete3_available=False, entrez_fallback=False):
    """Worker initializer: shared taxonomy cache and backend state (not inherited under spawn/forkserver)."""
    _taxonomy.cache = dict(cache_dict)
    _taxonomy.entrez_fallback = entrez_fallback
    if ete3_available and not _taxonomy.ete3_available:
        try:
            from ete3 import NCBITaxa

            _taxonomy.ncbi_taxa_cls = NCBITaxa
            _taxonomy.ete3_available = True
        except Exception as e:
            logger.warning(f"ete3 unavailable in worker: {e}")
