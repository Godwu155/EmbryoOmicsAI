import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite
import yaml

from embryoomics.config import InputError
from embryoomics.io import load_counts
from embryoomics.pipeline import inspect_input, run_analysis


@pytest.fixture
def sample(tmp_path):
    rng = np.random.default_rng(17)
    values = rng.poisson(1.1, size=(30, 25))
    values[:15, :5] += 5
    values[15:, 5:10] += 5
    values[:, 10] += 1
    folder = tmp_path / "tenx"
    folder.mkdir()
    matrix_file = folder / "GSMtest_matrix.mtx"
    mmwrite(matrix_file, sparse.csr_matrix(values.T))
    with matrix_file.open("rb") as source, gzip.open(folder / "GSMtest_matrix.mtx.gz", "wb") as target:
        target.write(source.read())
    matrix_file.unlink()
    (folder / "GSMtest_features.tsv").write_text("".join(f"id{i}\t{'mt-Nd1' if i == 0 else f'Gene{i}'}\tGene Expression\n" for i in range(25)), encoding="utf-8")
    (folder / "GSMtest_barcodes.tsv").write_text("".join(f"cell{i}\n" for i in range(30)), encoding="utf-8")
    config = {"dataset_accession": "SYNTHETIC", "sample_id": "test", "stage": "test_stage", "species": "mus_musculus", "genome_build": "unknown",
              "source_url": "https://example.invalid/synthetic", "raw_count_source": "X", "seed": 42,
              "qc": {"min_genes": 1, "max_genes": None, "max_pct_mt": None},
              "analysis": {"n_top_genes": 15, "n_neighbors": 5, "leiden_resolution": 0.5}}
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    return folder, config_path, config, values


def test_import_qc_and_metadata(sample, tmp_path):
    folder, config_path, config, values = sample
    info = inspect_input(folder, config_path)
    assert (info["cells"], info["genes"], info["retained"]) == (30, 25, 30)
    assert info["mt"]["matched_genes"] == 1
    assert info["qc_metrics"].iloc[0].total_counts == values[0].sum()
    metadata = tmp_path / "metadata.csv"
    pd.DataFrame({"cell_id": [f"cell{i}" for i in range(29)], "sample_id": "test", "stage": "test_stage"}).to_csv(metadata, index=False)
    with pytest.raises(InputError, match="METADATA_MISMATCH"):
        load_counts(folder, config, metadata)
    metadata = tmp_path / "metadata.csv"
    pd.DataFrame({"cell_id": ["cell0"] * 30}).to_csv(metadata, index=False)
    with pytest.raises(InputError, match="METADATA_DUPLICATE"):
        load_counts(folder, config, metadata)


def test_duplicate_ids_and_float_counts(sample, tmp_path):
    folder, _, config, values = sample
    (folder / "GSMtest_barcodes.tsv").write_text("cell0\n" * 30, encoding="utf-8")
    with pytest.raises(InputError, match="DUPLICATE_ID"):
        load_counts(folder, config)
    data = ad.AnnData(sparse.csr_matrix(values.astype(float) / 10))
    data.obs_names = [f"cell{i}" for i in range(30)]
    data.var_names = [f"id{i}" for i in range(25)]
    path = tmp_path / "normalized.h5ad"
    data.write_h5ad(path)
    with pytest.raises(InputError, match="COUNTS_NOT_RAW"):
        load_counts(path, config)


def test_no_mt_and_threshold_rejected(sample):
    folder, config_path, config, _ = sample
    feature_file = folder / "GSMtest_features.tsv"
    feature_file.write_text(feature_file.read_text(encoding="utf-8").replace("mt-Nd1", "Gene0"), encoding="utf-8")
    info = inspect_input(folder, config_path)
    assert info["mt"]["matched_genes"] == 0
    assert info["qc_metrics"].pct_counts_mt.isna().all()
    config["qc"]["max_pct_mt"] = 10
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    with pytest.raises(InputError, match="MT_UNAVAILABLE"):
        inspect_input(folder, config_path)


def test_end_to_end_preserves_counts(sample, tmp_path, monkeypatch):
    monkeypatch.delenv("EMBRYOOMICS_API_KEY", raising=False)
    folder, config_path, _, values = sample
    out = tmp_path / "run"
    result = run_analysis(config_path, folder, out)
    assert result["cells"] == 30
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert len(manifest["input_files"]) == 3
    assert all(len(entry["sha256"]) == 64 for entry in manifest["input_files"])
    assert manifest["seed"] == 42 and manifest["ai"]["enabled"] is False
    assert all((out / name).is_file() for name in manifest["files"])
    processed = ad.read_h5ad(out / "processed.h5ad")
    assert sparse.issparse(processed.layers["counts"])
    np.testing.assert_array_equal(processed.layers["counts"].toarray(), values)
    with pytest.raises(InputError, match="OUTPUT_EXISTS"):
        run_analysis(config_path, folder, out)
    with pytest.raises(InputError, match="OVERWRITE_UNSAFE"):
        run_analysis(config_path, folder, out, overwrite=True)
