"""Static HTML report and manifest, generated without an API key."""

from html import escape
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .i18n import columns_for, t, warning_for


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
    def table(frame, language):
        shown = frame.copy()
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
        result_text = (f"<p>{tr('umap_limit')}</p><figure><img src='figures/umap.png' alt='{tr('umap_figure')}'><figcaption>{tr('umap_figure')} {tr('umap_axes')}</figcaption></figure>"
                       + table(markers.head(100), language))
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
