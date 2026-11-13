"""Static HTML report and manifest, generated without an API key."""

from html import escape
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from .i18n import columns_for, t, warning_for


def save_marker_figures(data, out, markers):
    """Bound plotting to a small ranked gene panel; preserve the full sparse matrix."""
    figures = Path(out) / 'figures'
    figures.mkdir(exist_ok=True)
    # marker_table preserves Scanpy rank order within each cluster.
    candidates = markers.loc[markers.logfoldchange > 0] if not markers.empty else markers
    panel = candidates.groupby('cluster_id', sort=False).head(3).drop_duplicates('gene_id').head(30)
    panel = panel.loc[panel.gene_id.isin(data.var_names)].copy()
    files = ['marker_plot_genes.csv']
    panel.to_csv(Path(out) / files[0], index=False)
    details = {'backend': 'scanpy.pl.dotplot / scanpy.pl.matrixplot',
               'selection': 'First 3 positive-logfoldchange markers per cluster in markers.csv rank order; unique gene IDs; at most 30 genes',
               'gene_ids': panel.gene_id.tolist(), 'expression': 'X: log1p normalized expression; use_raw=False; no scaling',
               'detection': 'Fraction of cells with X > 0', 'files': files}
    if panel.empty:
        details['status'] = 'skipped_no_positive_markers'
        return details
    # Slice before Scanpy may densify; never materialize the full expression matrix.
    view = data[:, panel.gene_id.tolist()].copy()
    view.obs['cluster'] = pd.Categorical(view.obs['cluster'].astype(str))
    labels = {row.gene_id: f'{row.gene_symbol} ({row.gene_id})' for row in panel.itertuples()}
    view.var_names = [labels[gene] for gene in view.var_names]
    for name, plotter in [('marker_dotplot', sc.pl.dotplot), ('marker_heatmap', sc.pl.matrixplot)]:
        plot = plotter(view, var_names=view.var_names.tolist(), groupby='cluster', use_raw=False,
                       dendrogram=False, return_fig=True, show=False,
                       figsize=(max(7, len(panel) * 0.55), max(4, view.obs['cluster'].nunique() * 0.4 + 2)))
        try:
            plot.savefig(figures / f'{name}.png', dpi=140, bbox_inches='tight')
        finally:
            if plot.fig is not None:
                plt.close(plot.fig)
        files.append(f'figures/{name}.png')
    details['status'] = 'generated'
    return details


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


def save_doublet_figure(out, scores, simulated, threshold):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    try:
        for axis, values, title in zip(axes, [scores.doublet_score, simulated], ['Observed cells', 'Simulated doublets']):
            axis.hist(values, bins=40)
            axis.axvline(threshold, color='red', label=f'Threshold = {threshold:g}')
            axis.set(title=title, xlabel='Scrublet score', ylabel='Count')
            axis.legend()
        fig.tight_layout()
        fig.savefig(Path(out) / 'figures/doublet_scores.png', dpi=140)
    finally:
        plt.close(fig)


def write_report(out, manifest, qc_summary, markers, annotations):
    out = Path(out)
    def table(frame, language):
        shown = frame.copy()
        if 'filter_rule' in shown:
            shown['filter_rule'] = shown.filter_rule.map(lambda rule: t(language, 'qc_' + rule))
        if language == "zh-CN":
            replacements = {
                "test_method": {"Scanpy Wilcoxon on log1p normalized expression; cluster vs rest": "Scanpy Wilcoxon 检验：log1p 归一化表达，聚类与其余细胞比较"},
                "review_status": {"unreviewed": "未审核"},
                "candidate_cell_type": {"unknown": "未知"},
            }
            for column, values in replacements.items():
                if column in shown:
                    shown[column] = shown[column].replace(values)
        return shown.rename(columns=columns_for(language)).to_html(index=False, escape=True)

    def content(language):
        tr = lambda key: t(language, key)
        unknown_value = lambda value: tr("unknown") if value == "unknown" else escape(str(value))
        dimensions = manifest["input_dimensions"]
        data_text = (
            f"<p>{tr('dataset')}: {escape(str(manifest['dataset_accession']))}; "
            f"{tr('sample')}: {escape(str(manifest['sample_id']))}; "
            f"{tr('stage')}: {escape(str(manifest['stage']))}; "
            f"{tr('source')}: <a href='{escape(str(manifest['source_url']), quote=True)}'>{tr('geo_record')}</a>.</p>"
            f"<p>{tr('input_dimensions')}: {dimensions['cells']} {tr('cells')} × {dimensions['genes']} {tr('genes')}. "
            f"{tr('genome')}: {unknown_value(manifest['genome_build'])}. "
            f"{tr('license')}: {unknown_value(manifest['data_license'])}.</p>"
        )
        qc_text = table(qc_summary, language) + f"<figure><img src='figures/qc.png' alt='{tr('qc_figure')}'><figcaption>{tr('qc_figure')} {tr('qc_axes')}</figcaption></figure>"
        if manifest.get('qc', {}).get('filter_reasons'):
            qc_text += table(pd.DataFrame(manifest['qc']['filter_reasons']), language)
            qc_text += f"<p>{tr('qc_reasons_help')}</p>"
        doublets = manifest.get('doublets', {})
        if doublets.get('status') == 'completed':
            settings = doublets['settings']
            qc_text += (f"<h3>{tr('doublets_summary')}</h3><p>{tr('doublets_help')}</p>"
                        f"<p>{tr('doublets_mode')}: {tr('doublets_' + settings['mode'])}; "
                        f"{tr('doublets_rate')}: {settings['expected_doublet_rate']}; "
                        f"{tr('doublets_threshold')}: {settings['threshold']}; "
                        f"{tr('doublets_pcs')}: {settings['n_prin_comps']}; "
                        f"{tr('seed')}: {doublets['seed']}; "
                        f"{tr('doublets_upstream')}: {tr('doublets_upstream_' + settings['upstream_handling'])}; "
                        f"{tr('doublets_candidates')}: {doublets['candidate_doublets']}; "
                        f"{tr('qc_removed_cells')}: {doublets['removed_cells']}</p>"
                        f"<figure><img src='figures/doublet_scores.png' alt='{tr('doublets_summary')}'></figure>")
        result_text = (f"<p>{tr('umap_limit')}</p><figure><img src='figures/umap.png' alt='{tr('umap_figure')}'><figcaption>{tr('umap_figure')} {tr('umap_axes')}</figcaption></figure>"
                       + table(markers.head(100), language))
        for name in ('marker_dotplot', 'marker_heatmap'):
            if f'figures/{name}.png' in manifest['files']:
                result_text += f"<figure><img src='figures/{name}.png' alt='{tr(name)}'><figcaption>{tr(name)}</figcaption></figure>"
        if manifest.get('marker_plots', {}).get('status') == 'generated':
            result_text += f"<p>{tr('marker_plots_help')}</p>"
        review_text = f"<p>{tr('review_limit')}</p>" + table(annotations, language)
        sections = [("data", tr("data"), data_text), ("qc", tr("qc"), qc_text),
                    ("results", tr("report_results"), result_text), ("review", tr("report_review"), review_text)]
        nav = " ".join(f"<a href='#{language}-{identifier}'>{title}</a>" for identifier, title, _ in sections)
        warnings = "".join(f"<li>{escape(warning_for(language, warning))}</li>" for warning in manifest["warnings"])
        body = "".join(f"<section id='{language}-{identifier}'><h2>{title}</h2>{section}</section>" for identifier, title, section in sections)
        downloads = " ".join(f"<a href='{escape(name, quote=True)}'>{escape(name)}</a>" for name in manifest["files"] if name != "report.html")
        return (f"<div id='content-{language}' lang='{language}'{' hidden' if language == 'zh-CN' else ''}>"
                f"<nav>{nav}</nav><p>{tr('seed')}: {manifest['seed']}. {tr('ai_status')}</p>"
                f"<ul>{warnings}</ul>{body}<h2>{tr('downloads')}</h2>{downloads}</div>")

    html = ("<!doctype html><html lang='en'><meta charset='utf-8'><title>EmbryoOmics AI report / 报告</title>"
            "<style>body{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#24303b}"
            "nav a{margin-right:1rem}section{margin:3rem 0}img{max-width:100%}"
            "table{border-collapse:collapse;font-size:12px;display:block;overflow:auto}"
            "td,th{border:1px solid #ddd;padding:4px}th{background:#edf1f5}[hidden]{display:none!important}</style>"
            "<h1>EmbryoOmics AI v0.1</h1><label for='report-language'>Language / 语言</label> "
            "<select id='report-language' onchange='setReportLanguage(this.value)'>"
            "<option value='en'>English</option><option value='zh-CN'>简体中文</option></select>"
            + content("en") + content("zh-CN")
            + """<script>
function setReportLanguage(language) {
  if (language !== 'en' && language !== 'zh-CN') language = 'en';
  document.getElementById('content-en').hidden = language !== 'en';
  document.getElementById('content-zh-CN').hidden = language !== 'zh-CN';
  document.getElementById('report-language').value = language;
  document.documentElement.lang = language;
  try { localStorage.setItem('embryoomics-report-language', language); } catch (_) {}
}
try { setReportLanguage(localStorage.getItem('embryoomics-report-language') || 'en'); }
catch (_) { setReportLanguage('en'); }
</script></html>""")
    (out / "report.html").write_text(html, encoding="utf-8")
    (out / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
