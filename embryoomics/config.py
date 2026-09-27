"""Configuration contract shared by CLI and Streamlit."""

from pathlib import Path
import yaml


class InputError(ValueError):
    """A user-correctable input or configuration error."""


def load_config(path):
    with Path(path).open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise InputError("CONFIG: config.yaml must be a mapping")
    required = ["dataset_accession", "sample_id", "stage", "species", "genome_build", "source_url", "raw_count_source", "seed"]
    missing = [key for key in required if not config.get(key) and config.get(key) != 0]
    if missing:
        raise InputError(f"CONFIG: missing required fields: {', '.join(missing)}")
    if config["species"] != "mus_musculus":
        raise InputError("SPECIES: v0.1 requires species: mus_musculus")
    if not isinstance(config["seed"], int) or isinstance(config["seed"], bool):
        raise InputError("CONFIG: seed must be an integer")
    if config["raw_count_source"] != "X" and not str(config["raw_count_source"]).startswith("layers/") and config["raw_count_source"] != "raw/X":
        raise InputError("CONFIG: raw_count_source must be X, raw/X, or layers/<name>")
    qc = config.setdefault("qc", {})
    for key in ("min_genes", "max_genes"):
        value = qc.setdefault(key, None)
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise InputError(f"CONFIG: qc.{key} must be a nonnegative integer or null")
    if qc["min_genes"] is not None and qc["max_genes"] is not None and qc["min_genes"] > qc["max_genes"]:
        raise InputError("CONFIG: qc.min_genes exceeds qc.max_genes")
    value = qc.setdefault("max_pct_mt", None)
    if value is not None and (not isinstance(value, (int, float)) or not 0 <= value <= 100):
        raise InputError("CONFIG: qc.max_pct_mt must be 0–100 or null")
    analysis = config.setdefault("analysis", {})
    for key, default in (("n_top_genes", 2000), ("n_neighbors", 15)):
        value = analysis.setdefault(key, default)
        if not isinstance(value, int) or isinstance(value, bool) or value < 2:
            raise InputError(f"CONFIG: analysis.{key} must be an integer >= 2")
    value = analysis.setdefault("leiden_resolution", 0.5)
    if not isinstance(value, (int, float)) or value <= 0:
        raise InputError("CONFIG: analysis.leiden_resolution must be positive")
    return config
