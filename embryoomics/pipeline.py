"""One analysis service for CLI and Streamlit."""

from pathlib import Path
import platform
import shutil
import time
import importlib.metadata as metadata
import psutil

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from .ai import ai_status
from .config import InputError, load_config
from .evidence import load_annotations, marker_table
from .io import checksums, load_counts
from .qc import calculate_qc, preview_qc
from .reporting import save_figures, write_report


def inspect_input(input_path, config_path, metadata_path=None):
    config = load_config(config_path)
    data, files = load_counts(input_path, config, metadata_path)
    mt, warnings = calculate_qc(data)
    mask, summary = preview_qc(data, config["qc"])
    return {"config": config, "cells": data.n_obs, "genes": data.n_vars, "mt": mt,
            "warnings": warnings, "qc_metrics": data.obs[["sample_id", "stage", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy(),
            "preview": summary, "retained": int(mask.sum()), "files": [str(p) for p in files]}


def run_analysis(config_path, input_path, out, metadata_path=None, annotations_path=None, overwrite=False):
    start = time.perf_counter()
    config = load_config(config_path)
    out = Path(out)
    if out.exists() and not overwrite:
        raise InputError(f"OUTPUT_EXISTS: {out}; choose another path or --overwrite")
    data, files = load_counts(input_path, config, metadata_path)
    n_input = data.n_obs
    n_genes = data.n_vars
    mt, warnings = calculate_qc(data)
    metrics = data.obs[["sample_id", "stage", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    mask, qc_summary = preview_qc(data, config["qc"])
    if mask.sum() < 4:
        raise InputError(f"TOO_FEW_CELLS: {int(mask.sum())} cells remain; at least 4 are required for PCA/neighbors")
    data = data[mask].copy()
    if (data.obs["total_counts"] == 0).any():
        raise InputError("ZERO_COUNT_CELL: QC retained empty cells; set min_genes >= 1")
    if np.count_nonzero(np.asarray((data.layers["counts"] > 0).sum(axis=0)).ravel()) < 3:
        raise InputError("TOO_FEW_GENES: fewer than 3 expressed genes remain")
    # Keep untouched sparse integer counts in the output layer, while X is the working view.
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    n_top = min(config["analysis"]["n_top_genes"], data.n_vars)
    sc.pp.highly_variable_genes(data, n_top_genes=n_top, flavor="seurat")
    if int(data.var["highly_variable"].sum()) < 3:
        raise InputError("TOO_FEW_HVG: fewer than 3 highly variable genes; use a richer sample")
    model = data[:, data.var["highly_variable"]].copy()
    components = min(50, model.n_obs - 2, model.n_vars - 1)
    if components < 2:
        raise InputError("TOO_SMALL: PCA needs at least 2 components")
    sc.tl.pca(model, n_comps=components, svd_solver="arpack", random_state=config["seed"])
    neighbors = min(config["analysis"]["n_neighbors"], model.n_obs - 1)
    sc.pp.neighbors(model, n_neighbors=neighbors, random_state=config["seed"])
    sc.tl.leiden(model, resolution=config["analysis"]["leiden_resolution"], random_state=config["seed"], key_added="cluster", flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(model, random_state=config["seed"])
    data.obs["cluster"] = model.obs["cluster"].astype(str).to_numpy()
    data.obsm["X_pca"] = model.obsm["X_pca"]
    data.obsm["X_umap"] = model.obsm["X_umap"]
    data.uns["analysis_parameters"] = {"hvg_flavor": "seurat", "n_top_genes": n_top, "n_pcs": components, "n_neighbors": neighbors,
                                       "leiden_resolution": config["analysis"]["leiden_resolution"], "leiden_flavor": "igraph", "seed": config["seed"]}
    markers = marker_table(data)
    clusters = sorted(data.obs["cluster"].astype(str).unique())
    annotations = load_annotations(annotations_path, clusters)
    composition = data.obs.groupby(["cluster", "sample_id", "stage"], observed=True).size().rename("cells").reset_index()
    if "embryo_id" in data.obs:
        composition = data.obs.groupby(["cluster", "sample_id", "stage", "embryo_id"], observed=True).size().rename("cells").reset_index()
    warnings.append("Single-sample clustering and marker tests are exploratory; cell counts are not biological replicate counts.")
    if not config.get("genome_build") or config["genome_build"] == "unknown":
        warnings.append("Genome build is unknown; gene ID mapping is not verified.")
    if out.exists() and overwrite:
        shutil.rmtree(out)
    out.mkdir(parents=True)
    metrics.to_csv(out / "qc_metrics.csv", index_label="cell_id")
    qc_summary.to_csv(out / "qc_preview.csv", index=False)
    markers.to_csv(out / "markers.csv", index=False)
    annotations.to_csv(out / "annotations.csv", index=False)
    composition.to_csv(out / "cluster_composition.csv", index=False)
    (out / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")
    data.write_h5ad(out / "processed.h5ad", compression="gzip")
    save_figures(data, out, metrics)
    versions = {name: metadata.version(name) for name in ("scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib", "python-igraph", "leidenalg", "streamlit", "psutil")}
    memory = psutil.Process().memory_info()
    peak_memory = getattr(memory, "peak_wset", None)
    manifest = {"dataset_accession": config["dataset_accession"], "sample_id": config["sample_id"], "stage": config["stage"],
                "species": config["species"], "genome_build": config["genome_build"], "source_url": config["source_url"],
                "raw_count_source": config["raw_count_source"], "input_files": checksums(files),
                "input_dimensions": {"cells": n_input, "genes": n_genes}, "qc": {"input_cells": n_input, "retained_cells": data.n_obs, "mt_detection": mt, "preview_by_sample": qc_summary.to_dict(orient="records")},
                "parameters": config, "seed": config["seed"], "versions": {"python": platform.python_version(), **versions},
                "warnings": warnings, "ai": ai_status(), "elapsed_seconds": round(time.perf_counter() - start, 2),
                "peak_working_set_bytes": peak_memory,
                "files": ["report.html", "run_manifest.json", "qc_metrics.csv", "qc_preview.csv", "markers.csv", "annotations.csv", "cluster_composition.csv", "processed.h5ad", "config.yaml", "figures/qc.png", "figures/umap.png"]}
    write_report(out, manifest, qc_summary, markers, annotations)
    return {"output": str(out.resolve()), "cells": data.n_obs, "genes": data.n_vars, "clusters": len(clusters), "elapsed_seconds": manifest["elapsed_seconds"]}
