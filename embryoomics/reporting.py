"""Static HTML report and manifest, generated without an API key."""

from html import escape
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def save_figures(data, out, input_qc_metrics):
    figures = Path(out) / "figures"
    figures.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    for axis, key, title in zip(axes, ["total_counts", "n_genes_by_counts", "pct_counts_mt"], ["Total counts", "Detected genes", "Mitochondrial %"]):
        values = input_qc_metrics[key].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if len(values):
            axis.hist(values, bins=min(40, max(5, len(values) // 5)))
        axis.set_title(title)
    fig.tight_layout()
    fig.savefig(figures / "qc.png", dpi=140)
    plt.close(fig)
    if "X_umap" in data.obsm:
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
        coords = data.obsm["X_umap"]
        for axis, key in zip(axes, ["cluster", "stage", "sample_id"]):
            categories = data.obs[key].astype(str)
            for value in sorted(categories.unique()):
                selected = (categories == value).to_numpy()
                axis.scatter(coords[selected, 0], coords[selected, 1], s=5, label=value, rasterized=True)
            axis.set_title(key)
            axis.legend(fontsize=6, loc="best")
            axis.set_xticks([])
            axis.set_yticks([])
        fig.tight_layout()
        fig.savefig(figures / "umap.png", dpi=140)
        plt.close(fig)


def write_report(out, manifest, qc_summary, markers, annotations):
    out = Path(out)
    sections = [
        ("data", "Data", f"<p>Dataset: {escape(manifest['dataset_accession'])}; sample: {escape(manifest['sample_id'])}; stage: {escape(manifest['stage'])}; source: <a href='{escape(manifest['source_url'], quote=True)}'>GEO record</a>.</p><p>Input dimensions: {manifest['input_dimensions']['cells']} cells × {manifest['input_dimensions']['genes']} genes. Genome: {escape(manifest['genome_build'])}.</p>"),
        ("qc", "QC", qc_summary.to_html(index=False, escape=True) + "<img src='figures/qc.png' alt='QC distributions'>"),
        ("results", "Results and marker evidence", "<p>UMAP shows similarity in this analysis; it is not a developmental trajectory.</p><img src='figures/umap.png' alt='UMAP by cluster, stage and sample'>" + markers.head(100).to_html(index=False, escape=True)),
        ("report", "Review and limitations", "<p>Cell type labels require human review. Cluster markers are exploratory and cells are not biological replicates. No lineage, causal, or cross-modality claim is made.</p>" + annotations.to_html(index=False, escape=True)),
    ]
    nav = " ".join(f"<a href='#{identifier}'>{title}</a>" for identifier, title, _ in sections)
    warnings = "".join(f"<li>{escape(w)}</li>" for w in manifest["warnings"])
    body = "".join(f"<section id='{identifier}'><h2>{title}</h2>{content}</section>" for identifier, title, content in sections)
    downloads = " ".join(f"<a href='{escape(name, quote=True)}'>{escape(name)}</a>" for name in manifest["files"] if name != "report.html")
    html = f"""<!doctype html><html lang='en'><meta charset='utf-8'><title>EmbryoOmics AI report</title><style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#24303b}}nav a{{margin-right:1rem}}section{{margin:3rem 0}}img{{max-width:100%}}table{{border-collapse:collapse;font-size:12px;display:block;overflow:auto}}td,th{{border:1px solid #ddd;padding:4px}}th{{background:#edf1f5}}</style><h1>EmbryoOmics AI v0.1</h1><nav>{nav}</nav><p>Analysis seed: {manifest['seed']}. AI interpretation: disabled.</p><ul>{warnings}</ul>{body}<h2>Downloads</h2>{downloads}</html>"""
    (out / "report.html").write_text(html, encoding="utf-8")
    (out / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
