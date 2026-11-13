"""Configuration contract shared by CLI and Streamlit."""

from pathlib import Path
from copy import deepcopy
import math
import yaml


class InputError(ValueError):
    """A user-correctable input or configuration error."""


def load_config(path):
    with Path(path).open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    return validate_config(config)


def validate_config(config):
    config = deepcopy(config)
    if not isinstance(config, dict):
        raise InputError("CONFIG: config.yaml must be a mapping")
    required = ["dataset_accession", "sample_id", "stage", "species", "genome_build", "source_url", "raw_count_source", "seed"]
    missing = [key for key in required if not config.get(key) and config.get(key) != 0]
    if missing:
        raise InputError(f"CONFIG: missing required fields: {', '.join(missing)}")
    if config["species"] not in ("mus_musculus", "homo_sapiens"):
        raise InputError("SPECIES: supported values are mus_musculus and homo_sapiens")
    if config.get("count_matrix_layout") not in (None, "genes_by_cells_tsv"):
        raise InputError("CONFIG: count_matrix_layout must be genes_by_cells_tsv for a tabular count matrix")
    config.setdefault("data_license", "unknown")
    mt_ids = config.setdefault('mitochondrial_gene_ids', None)
    if mt_ids is not None and (not isinstance(mt_ids, list) or not mt_ids or
            any(not isinstance(gene, str) or not gene.strip() for gene in mt_ids) or len(set(mt_ids)) != len(mt_ids)):
        raise InputError('CONFIG: mitochondrial_gene_ids must be a nonempty list of unique gene IDs or null')
    if not isinstance(config["seed"], int) or isinstance(config["seed"], bool):
        raise InputError("CONFIG: seed must be an integer")
    if config["raw_count_source"] != "X" and not str(config["raw_count_source"]).startswith("layers/") and config["raw_count_source"] != "raw/X":
        raise InputError("CONFIG: raw_count_source must be X, raw/X, or layers/<name>")
    qc = config.setdefault("qc", {})
    if not isinstance(qc, dict):
        raise InputError("CONFIG: qc must be a mapping")
    for key in ("min_genes", "max_genes"):
        value = qc.setdefault(key, None)
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise InputError(f"CONFIG: qc.{key} must be a nonnegative integer or null")
    if qc["min_genes"] is not None and qc["max_genes"] is not None and qc["min_genes"] > qc["max_genes"]:
        raise InputError("CONFIG: qc.min_genes exceeds qc.max_genes")
    value = qc.setdefault("max_pct_mt", None)
    if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 100):
        raise InputError("CONFIG: qc.max_pct_mt must be 0–100 or null")
    analysis = config.setdefault("analysis", {})
    if not isinstance(analysis, dict):
        raise InputError("CONFIG: analysis must be a mapping")
    for key, default in (("n_top_genes", 2000), ("n_neighbors", 15)):
        value = analysis.setdefault(key, default)
        if not isinstance(value, int) or isinstance(value, bool) or value < 2:
            raise InputError(f"CONFIG: analysis.{key} must be an integer >= 2")
    value = analysis.setdefault("leiden_resolution", 0.5)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise InputError("CONFIG: analysis.leiden_resolution must be positive")
    doublets = config.setdefault('doublets', {})
    if not isinstance(doublets, dict):
        raise InputError('CONFIG: doublets must be a mapping')
    if doublets.setdefault('mode', 'off') not in ('off', 'score', 'filter'):
        raise InputError('CONFIG: doublets.mode must be off, score or filter')
    for key in ('expected_doublet_rate', 'threshold'):
        value = doublets.setdefault(key, None)
        if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value < 1):
            raise InputError(f'CONFIG: doublets.{key} must be between 0 and 1, exclusively')
        if doublets['mode'] != 'off' and value is None:
            raise InputError(f'CONFIG: set doublets.{key} explicitly when enabled')
    value = doublets.setdefault('n_prin_comps', 30)
    if not isinstance(value, int) or isinstance(value, bool) or value < 2:
        raise InputError('CONFIG: doublets.n_prin_comps must be an integer >= 2')
    if doublets.setdefault('upstream_handling', 'unknown') not in ('unknown', 'not_removed', 'removed'):
        raise InputError('CONFIG: doublets.upstream_handling must be unknown, not_removed or removed')
    return config
