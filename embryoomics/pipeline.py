"""One analysis service for CLI and Streamlit."""

from pathlib import Path
import platform
import shutil
import time
import importlib.metadata as metadata
import zipfile
import psutil

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from .ai import ai_status
from .config import InputError, load_config, validate_config
from .evidence import load_annotations, marker_table
from .io import checksums, load_counts, review_identity
from .qc import calculate_qc, preview_qc_metrics
from .doublets import score_doublets, apply_doublet_decision
from .reporting import save_figures, save_marker_figures, save_doublet_figure, write_report


def inspect_input(input_path, config_path, metadata_path=None, *, apply_config_filters=True):
    config = load_config(config_path)
    identity = review_identity(input_path, config, metadata_path)
    data, files = load_counts(input_path, config, metadata_path)
    mt, warnings = calculate_qc(data)
    thresholds = config['qc'] if apply_config_filters else {'min_genes': None, 'max_genes': None, 'max_pct_mt': None}
    mask, summary, reasons = preview_qc_metrics(data.obs, thresholds)
    if identity != review_identity(input_path, load_config(config_path), metadata_path):
        raise InputError('INPUT_CHANGED: input or configuration changed during validation; validate again')
    return {"config": config, "cells": data.n_obs, "genes": data.n_vars, "mt": mt,
            "warnings": warnings, "qc_metrics": data.obs[["sample_id", "stage", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy(),
            "preview": summary, "filter_reasons": reasons, "input_identity": identity,
            "retained": int(mask.sum()), "files": [str(p) for p in files]}


def inspect_doublets(input_path, config_path, metadata_path=None, *, qc_override=None, doublet_override=None):
    source = load_config(config_path)
    identity = review_identity(input_path, source, metadata_path)
    config = validate_config({**source, 'qc': qc_override if qc_override is not None else source['qc'],
                              'doublets': doublet_override if doublet_override is not None else source['doublets']})
    if config['doublets']['mode'] == 'off':
        raise InputError('SCRUBLET_OFF: enable scoring before requesting a preview')
    data, files = load_counts(input_path, config, metadata_path)
    calculate_qc(data)
    mask, _, _ = preview_qc_metrics(data.obs, config['qc'])
    scores, simulated, summary = score_doublets(data[mask], config['doublets'], config['seed'])
    final, audit = apply_doublet_decision(data.obs, mask, scores, config['doublets'])
    if identity != review_identity(input_path, load_config(config_path), metadata_path):
        raise InputError('QC_STALE: input or configuration changed while scoring; validate again')
    return {'scores': scores, 'simulated': simulated, 'summary': summary, 'audit': audit,
            'retained': int(final.sum()), 'input_identity': identity, 'config': config}


def run_analysis(config_path, input_path, out, metadata_path=None, annotations_path=None, overwrite=False,
                 *, qc_override=None, expected_input_identity=None, doublet_override=None,
                 confirm_doublets=False, expected_doublet_digest=None):
    start = time.perf_counter()
    config = load_config(config_path)
    identity = review_identity(input_path, config, metadata_path)
    if expected_input_identity is not None and identity != expected_input_identity:
        raise InputError('QC_STALE: input or configuration changed; validate and review QC again')
    if qc_override is not None:
        config = validate_config({**config, 'qc': qc_override})
    if doublet_override is not None:
        config = validate_config({**config, 'doublets': doublet_override})
    if config['doublets']['mode'] == 'filter' and not confirm_doublets:
        raise InputError('DOUBLET_CONFIRMATION: review doublet scores and explicitly confirm filtering')
    out = Path(out)
    if out.exists() and not overwrite:
        raise InputError(f"OUTPUT_EXISTS: {out}; choose another path or --overwrite")
    if out.exists() and overwrite:
        resolved = out.resolve()
        workspace = Path.cwd().resolve()
        if not out.is_dir() or workspace not in resolved.parents or not (resolved / "run_manifest.json").is_file():
            raise InputError("OVERWRITE_UNSAFE: only an existing EmbryoOmics run directory inside the current workspace may be overwritten")
    data, files = load_counts(input_path, config, metadata_path)
    if identity != review_identity(input_path, load_config(config_path), metadata_path):
        raise InputError('QC_STALE: input or configuration changed while loading; validate again')
    input_checksums = checksums(files)
    if identity != review_identity(input_path, load_config(config_path), metadata_path):
        raise InputError('QC_STALE: input or configuration changed during checksum calculation; validate again')
    n_input = data.n_obs
    n_genes = data.n_vars
    mt, warnings = calculate_qc(data)
    metrics = data.obs[["sample_id", "stage", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    mask, qc_summary, filter_reasons = preview_qc_metrics(data.obs, config["qc"])
    doublet_summary = {'settings': config['doublets'], 'status': 'disabled'}
    doublet_audit = None
    if config['doublets']['mode'] != 'off':
        scores, simulated, doublet_summary = score_doublets(data[mask], config['doublets'], config['seed'])
        if expected_doublet_digest is not None and expected_doublet_digest != doublet_summary['score_digest']:
            raise InputError('DOUBLET_STALE: scoring differs from the reviewed preview; preview again')
        ordinary_retained = int(mask.sum())
        mask, doublet_audit = apply_doublet_decision(data.obs, mask, scores, config['doublets'])
        doublet_summary.update(status='completed', filtering_confirmed=bool(confirm_doublets) if config['doublets']['mode'] == 'filter' else False)
        qc_summary['retained_cells'] = int(mask.sum())
        qc_summary['removed_cells'] = qc_summary.input_cells - qc_summary.retained_cells
        doublet_reason = {'filter_rule': 'doublets', 'threshold': config['doublets']['threshold'],
            'enabled': config['doublets']['mode'] == 'filter', 'failed_cells': int(scores.predicted_doublet.sum()),
            'newly_removed_cells': ordinary_retained - int(mask.sum()), 'retained_after_rule': int(mask.sum())}
        filter_reasons = pd.DataFrame([*filter_reasons.to_dict(orient='records'), doublet_reason])
        warnings.extend(doublet_summary['warnings'])
        warnings.append('Scrublet candidates are exploratory; upstream doublet handling may be unknown and labels are not confirmed doublets.')
    if mask.sum() < 4:
        raise InputError(f"TOO_FEW_CELLS: {int(mask.sum())} cells remain; at least 4 are required for PCA/neighbors")
    data = data[mask].copy()
    if doublet_audit is not None:
        for column in scores:
            data.obs[column] = scores.loc[data.obs_names, column].to_numpy()
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
        # The absolute target was resolved and checked against the workspace above.
        shutil.rmtree(out)
    out.mkdir(parents=True)
    metrics.to_csv(out / "qc_metrics.csv", index_label="cell_id")
    if doublet_audit is not None:
        doublet_audit.to_csv(out / 'doublet_scores.csv')
        pd.DataFrame({'doublet_score': simulated}).to_csv(out / 'doublet_simulated_scores.csv', index=False)
    qc_summary.to_csv(out / "qc_preview.csv", index=False)
    filter_reasons.to_csv(out / 'qc_filter_reasons.csv', index=False)
    markers.to_csv(out / "markers.csv", index=False)
    annotations.to_csv(out / "annotations.csv", index=False)
    composition.to_csv(out / "cluster_composition.csv", index=False)
    (out / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")
    data.write_h5ad(out / "processed.h5ad", compression="gzip")
    save_figures(data, out, metrics)
    if doublet_audit is not None:
        save_doublet_figure(out, scores, simulated, config['doublets']['threshold'])
    marker_plots = save_marker_figures(data, out, markers)
    versions = {name: metadata.version(name) for name in ("scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib", "python-igraph", "leidenalg", "streamlit", "psutil", "scikit-learn")}
    memory = psutil.Process().memory_info()
    peak_memory = getattr(memory, "peak_wset", None)
    manifest = {"dataset_accession": config["dataset_accession"], "sample_id": config["sample_id"], "stage": config["stage"],
                "species": config["species"], "genome_build": config["genome_build"], "source_url": config["source_url"],
                "data_license": config["data_license"],
                "raw_count_source": config["raw_count_source"], "input_files": input_checksums,
                "input_dimensions": {"cells": n_input, "genes": n_genes}, "qc": {"input_cells": n_input, "retained_cells": data.n_obs, "mt_detection": mt, "preview_by_sample": qc_summary.to_dict(orient="records"),
                    "filter_reasons": filter_reasons.astype(object).where(pd.notna(filter_reasons), None).to_dict(orient='records')},
                "parameters": config, "seed": config["seed"], "versions": {"python": platform.python_version(), **versions},
                "warnings": warnings, "ai": ai_status(), "elapsed_seconds": round(time.perf_counter() - start, 2),
                "peak_working_set_bytes": peak_memory,
                "files": ["report.html", "run_manifest.json", "qc_metrics.csv", "qc_preview.csv", "qc_filter_reasons.csv", "markers.csv", "annotations.csv", "cluster_composition.csv", "processed.h5ad", "config.yaml", "figures/qc.png", "figures/umap.png"]}
    manifest['marker_plots'] = marker_plots
    manifest['doublets'] = doublet_summary
    if doublet_audit is not None:
        manifest['files'].extend(['doublet_scores.csv', 'doublet_simulated_scores.csv', 'figures/doublet_scores.png'])
    manifest['files'].extend(marker_plots['files'])
    write_report(out, manifest, qc_summary, markers, annotations)
    # Keep relative image and download links usable after downloading the report.
    with zipfile.ZipFile(out / 'report_bundle.zip', 'w', compression=zipfile.ZIP_DEFLATED) as bundle:
        for name in manifest['files']:
            bundle.write(out / name, arcname=name)
    return {"output": str(out.resolve()), "cells": data.n_obs, "genes": data.n_vars, "clusters": len(clusters), "elapsed_seconds": manifest["elapsed_seconds"]}
