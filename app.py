"""Local Streamlit review UI using the same analysis service as CLI."""

from pathlib import Path

import pandas as pd
import numpy as np
import streamlit as st
import yaml

from embryoomics.config import InputError, load_config
from embryoomics.i18n import LANGUAGES, columns_for, t, warning_for
from embryoomics.pipeline import inspect_input, run_analysis


st.set_page_config(page_title="EmbryoOmics AI v0.1", layout="wide")
language = LANGUAGES[st.sidebar.selectbox("Language / 语言", list(LANGUAGES), key="language")]
st.title("EmbryoOmics AI v0.1")
st.caption(t(language, "subtitle"))
pages = st.tabs([t(language, key) for key in ("data", "qc", "results", "report")])

with st.sidebar:
    input_path = st.text_input(t(language, "input_path"), "data/GSM8559287", key="input_path")
    config_path = st.text_input(t(language, "config_path"), "examples/GSM8559287.yaml", key="config_path")
    metadata_path = st.text_input(t(language, "metadata_path"), "", key="metadata_path")
    out_path = st.text_input(t(language, "out_path"), "results/run_001", key="out_path")

try:
    config = load_config(config_path)
except (InputError, FileNotFoundError, yaml.YAMLError) as exc:
    st.error(f"{t(language, 'configuration')}: {exc}")
    st.stop()

with pages[0]:
    fields = ("dataset_accession", "sample_id", "stage", "species", "genome_build", "source_url", "data_license", "raw_count_source")
    labels = columns_for(language)
    st.dataframe(pd.DataFrame([{t(language, "field"): labels.get(key, key),
                                t(language, "value"): t(language, "unknown") if config[key] == "unknown" else config[key]} for key in fields]), hide_index=True)
    if st.button(t(language, "validate"), key="validate"):
        try:
            st.session_state["inspection"] = inspect_input(input_path, config_path, metadata_path or None)
            st.session_state["inspected_paths"] = (input_path, config_path, metadata_path)
        except (InputError, FileNotFoundError, ValueError) as exc:
            st.error(str(exc))
    if "inspection" in st.session_state:
        info = st.session_state["inspection"]
        st.success(f"{info['cells']} {t(language, 'cells')} × {info['genes']} {t(language, 'genes')}")
        mt_method = info["mt"]["method"]
        mt_method = t(language, "mt_symbol_method" if mt_method.startswith("mouse gene symbols") else "mt_no_match")
        st.write(t(language, "mt_detection"), {t(language, "mt_method"): mt_method,
                                               t(language, "mt_matched"): info["mt"]["matched_genes"]})
        for warning in info["warnings"]:
            st.warning(warning_for(language, warning))

with pages[1]:
    st.write(t(language, "qc_intro"))
    if "inspection" in st.session_state and st.session_state.get("inspected_paths") == (input_path, config_path, metadata_path):
        info = st.session_state["inspection"]
        st.dataframe(info["preview"].rename(columns=labels), hide_index=True)
        metrics = info["qc_metrics"]
        for key in ("total_counts", "n_genes_by_counts", "pct_counts_mt"):
            if metrics[key].notna().any():
                st.subheader(labels.get(key, key))
                counts, edges = np.histogram(metrics[key].dropna().to_numpy(dtype=float), bins=30)
                histogram = pd.DataFrame({"bin_midpoint": (edges[:-1] + edges[1:]) / 2, t(language, "chart_count"): counts}).set_index("bin_midpoint")
                st.bar_chart(histogram)
        st.write(t(language, "qc_thresholds"))
        st.dataframe(pd.DataFrame([{t(language, "field"): labels.get(key, key),
                                    t(language, "value"): value if value is not None else t(language, "unset")}
                                   for key, value in config["qc"].items()]), hide_index=True)
        st.caption(t(language, "qc_edit"))
        confirmed = st.checkbox(t(language, "qc_confirm"), key="qc_confirm")
        if st.button(t(language, "run"), disabled=not confirmed, key="run"):
            try:
                result = run_analysis(config_path, input_path, out_path, metadata_path or None)
                st.session_state["result"] = result
                st.success(f"{t(language, 'completed')}: {result['output']}")
            except (InputError, FileNotFoundError, ValueError) as exc:
                st.error(str(exc))
    else:
        st.info(t(language, "validate_first"))

result_dir = Path(st.session_state.get("result", {}).get("output", out_path))
with pages[2]:
    if (result_dir / "markers.csv").is_file():
        st.image(str(result_dir / "figures" / "umap.png"))
        st.caption(f"{t(language, 'umap_axes')} {t(language, 'umap_limit')}")
        markers = pd.read_csv(result_dir / "markers.csv")
        shown_markers = markers.copy()
        if language == "zh-CN" and "test_method" in shown_markers:
            shown_markers["test_method"] = shown_markers["test_method"].replace({
                "Scanpy Wilcoxon on log1p normalized expression; cluster vs rest": "Scanpy Wilcoxon 检验：log1p 归一化表达，聚类与其余细胞比较"})
        st.dataframe(shown_markers.rename(columns=labels))
        st.caption(t(language, "annotation_hint"))
        current = pd.read_csv(result_dir / "annotations.csv", dtype=str, keep_default_na=False)
        edited = st.data_editor(current, num_rows="fixed", hide_index=True,
                                column_config={key: st.column_config.TextColumn(label) for key, label in labels.items() if key in current})
        st.download_button(t(language, "annotation_download"), edited.to_csv(index=False).encode("utf-8"), "annotations_reviewed.csv", "text/csv", key="annotations_download")
        st.caption(t(language, "annotation_rerun"))
    else:
        st.info(t(language, "run_first"))

with pages[3]:
    if (result_dir / "report.html").is_file():
        st.write(f"{t(language, 'report_location')}: {result_dir / 'report.html'}")
        for filename in ("report.html", "run_manifest.json", "qc_metrics.csv", "markers.csv", "annotations.csv", "processed.h5ad"):
            st.download_button(filename, (result_dir / filename).read_bytes(), filename, key=filename)
        st.caption(t(language, "ai_disabled"))
    else:
        st.info(t(language, "report_first"))
