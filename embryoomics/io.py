"""Sparse counts loading, strict validation, and input provenance."""

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from .config import InputError


def input_files(path, metadata_path=None):
    path = Path(path)
    if not path.exists():
        raise InputError(f"INPUT_MISSING: {path} does not exist")
    if path.is_dir():
        matches = {}
        for suffix in ("matrix.mtx", "features.tsv", "barcodes.tsv"):
            found = sorted(p for p in path.iterdir() if p.name.endswith((suffix, suffix + ".gz")))
            if len(found) != 1:
                raise InputError(f"TENX_FILES: expected exactly one *{suffix}[.gz] in {path}, found {len(found)}")
            matches[suffix] = found[0]
        prefixes = [matches[s].name.split(s)[0] for s in matches]
        if len(set(prefixes)) != 1:
            raise InputError("TENX_FILES: matrix/features/barcodes prefixes differ; select one sample only")
        files = list(matches.values())
    elif path.suffix.lower() == ".h5ad":
        files = [path]
    else:
        raise InputError("FORMAT: provide a 10x three-file folder or .h5ad; CSV import is not enabled")
    if metadata_path:
        meta = Path(metadata_path)
        if not meta.is_file():
            raise InputError(f"METADATA_MISSING: {meta} does not exist")
        files.append(meta)
    return files


def checksums(files):
    result = []
    for path in files:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        result.append({"file": str(Path(path).resolve()), "bytes": Path(path).stat().st_size, "sha256": digest.hexdigest()})
    return result


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", header=None, dtype=str, compression="infer", keep_default_na=False)


def load_counts(path, config, metadata_path=None):
    files = input_files(path, metadata_path)
    path = Path(path)
    if path.is_dir():
        if config["raw_count_source"] != "X":
            raise InputError("COUNT_SOURCE: 10x folder requires raw_count_source: X")
        selected = {suffix: next(p for p in files if p.name.endswith((suffix, suffix + ".gz"))) for suffix in ("matrix.mtx", "features.tsv", "barcodes.tsv")}
        matrix = mmread(str(selected["matrix.mtx"]))
        if not sparse.issparse(matrix):
            raise InputError("SPARSE_REQUIRED: Matrix Market input must be sparse")
        matrix = matrix.tocsr().T.tocsr()
        features = _read_tsv(selected["features.tsv"])
        barcodes = _read_tsv(selected["barcodes.tsv"])
        if features.shape[1] < 2 or barcodes.shape[1] != 1:
            raise InputError("TENX_COLUMNS: features needs ID and symbol; barcodes needs one column")
        if matrix.shape != (len(barcodes), len(features)):
            raise InputError(f"DIMENSIONS: matrix cells×genes {matrix.shape} differs from barcodes/features {(len(barcodes), len(features))}")
        data = ad.AnnData(matrix)
        data.obs_names = pd.Index(barcodes.iloc[:, 0].astype(str), name="cell_id")
        data.var_names = pd.Index(features.iloc[:, 0].astype(str), name="gene_id")
        data.var["gene_symbol"] = features.iloc[:, 1].to_numpy(dtype=str)
        if features.shape[1] >= 3:
            data.var["feature_type"] = features.iloc[:, 2].to_numpy(dtype=str)
        if "feature_type" in data.var and not (data.var["feature_type"] == "Gene Expression").all():
            raise InputError("FEATURE_TYPE: folder includes non-gene-expression features; prepare an RNA-only triplet")
    else:
        data = ad.read_h5ad(path)
        source = config["raw_count_source"]
        if source == "X":
            matrix = data.X
        elif source == "raw/X":
            if data.raw is None:
                raise InputError("COUNT_SOURCE: raw/X is absent")
            raw = data.raw.to_adata()
            data = raw
            matrix = raw.X
        else:
            layer = source.split("/", 1)[1]
            if layer not in data.layers:
                raise InputError(f"COUNT_SOURCE: layer {layer!r} is absent")
            matrix = data.layers[layer]
        data.X = matrix
        if "gene_symbol" not in data.var:
            data.var["gene_symbol"] = data.var_names.astype(str)
    if data.n_obs == 0 or data.n_vars == 0:
        raise InputError("EMPTY: matrix needs at least one cell and gene")
    for axis, ids in (("cell", data.obs_names), ("gene", data.var_names)):
        values = pd.Index(ids.astype(str))
        duplicates = values[values.duplicated()].unique().tolist()
        if duplicates or (values == "").any():
            raise InputError(f"DUPLICATE_ID: {axis} IDs must be nonempty and unique; examples: {duplicates[:5]}")
    if not sparse.issparse(data.X):
        raise InputError("SPARSE_REQUIRED: counts must be stored sparsely")
    counts = data.X.tocsr(copy=True)
    counts.sum_duplicates()
    if not np.issubdtype(counts.dtype, np.number) or not np.isfinite(counts.data).all() or (counts.data < 0).any() or not np.equal(counts.data, np.floor(counts.data)).all():
        raise InputError("COUNTS_NOT_RAW: counts must be finite nonnegative integers; choose the raw counts layer")
    counts.eliminate_zeros()
    data.X = counts.astype(np.int64)
    if metadata_path:
        metadata = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
        if "cell_id" not in metadata:
            raise InputError("METADATA: metadata.csv needs cell_id")
        duplicate = metadata.loc[metadata.cell_id.duplicated(), "cell_id"].head().tolist()
        if duplicate:
            raise InputError(f"METADATA_DUPLICATE: duplicate cell_id examples: {duplicate}")
        given = set(metadata.cell_id)
        expected = set(data.obs_names)
        if given != expected:
            raise InputError(f"METADATA_MISMATCH: {len(expected-given)} missing cells, {len(given-expected)} extra; examples: {list(expected-given)[:3]}, {list(given-expected)[:3]}")
        metadata = metadata.set_index("cell_id").loc[data.obs_names]
        for key in ("sample_id", "stage"):
            if key not in metadata:
                metadata[key] = config[key]
            elif (metadata[key] == "").any():
                raise InputError(f"METADATA: empty {key} is not allowed")
        if (metadata["sample_id"] != config["sample_id"]).any() or (metadata["stage"] != config["stage"]).any():
            raise InputError("METADATA_CONFLICT: sample_id/stage differs from single-sample config")
        for key in ("sample_id", "stage", "embryo_id", "batch"):
            if key in metadata:
                data.obs[key] = metadata[key].to_numpy(dtype=str)
    else:
        data.obs["sample_id"] = config["sample_id"]
        data.obs["stage"] = config["stage"]
    data.layers["counts"] = data.X.copy()
    return data, files
