"""Local Streamlit review UI using the same analysis service as CLI."""

from pathlib import Path

import pandas as pd
import streamlit as st
import yaml

from embryoomics.config import InputError, load_config
from embryoomics.pipeline import inspect_input, run_analysis


st.set_page_config(page_title="EmbryoOmics AI v0.1", layout="wide")
st.title("EmbryoOmics AI v0.1")
st.caption("Local mouse embryo scRNA-seq review. Cell types need human confirmation; UMAP is not a trajectory.")
pages = st.tabs(["Data", "QC", "Results", "Report"])

with st.sidebar:
    input_path = st.text_input("10x folder or .h5ad path", "data/GSM8559287")
    config_path = st.text_input("Config YAML path", "examples/GSM8559287.yaml")
    metadata_path = st.text_input("Optional metadata CSV path", "")
    out_path = st.text_input("New results folder", "results/run_001")

try:
    config = load_config(config_path)
except (InputError, FileNotFoundError, yaml.YAMLError) as exc:
    st.error(f"Configuration: {exc}")
    st.stop()

with pages[0]:
    st.write({key: config[key] for key in ("dataset_accession", "sample_id", "stage", "species", "genome_build", "source_url", "raw_count_source")})
    if st.button("Validate input"):
        try:
            st.session_state["inspection"] = inspect_input(input_path, config_path, metadata_path or None)
            st.session_state["inspected_paths"] = (input_path, config_path, metadata_path)
        except (InputError, FileNotFoundError, ValueError) as exc:
            st.error(str(exc))
    if "inspection" in st.session_state:
        info = st.session_state["inspection"]
        st.success(f"{info['cells']} cells × {info['genes']} genes")
        st.write("Mitochondrial detection", info["mt"])
        for warning in info["warnings"]:
            st.warning(warning)

with pages[1]:
    st.write("Review QC distributions and choose thresholds. A run requires explicit confirmation.")
    if "inspection" in st.session_state and st.session_state.get("inspected_paths") == (input_path, config_path, metadata_path):
        info = st.session_state["inspection"]
        st.dataframe(info["preview"], hide_index=True)
        metrics = info["qc_metrics"]
        for key in ("total_counts", "n_genes_by_counts", "pct_counts_mt"):
            if metrics[key].notna().any():
                st.subheader(key)
                st.bar_chart(metrics[key].value_counts(bins=30, sort=False))
        st.write("Configured thresholds", config["qc"])
        st.caption("Edit the YAML, then validate again to update the preview.")
        confirmed = st.checkbox("I reviewed the QC preview and confirm these thresholds")
        if st.button("Run analysis", disabled=not confirmed):
            try:
                result = run_analysis(config_path, input_path, out_path, metadata_path or None)
                st.session_state["result"] = result
                st.success(f"Completed: {result['output']}")
            except (InputError, FileNotFoundError, ValueError) as exc:
                st.error(str(exc))
    else:
        st.info("Validate input on the Data tab first.")

result_dir = Path(st.session_state.get("result", {}).get("output", out_path))
with pages[2]:
    if (result_dir / "markers.csv").is_file():
        st.image(str(result_dir / "figures" / "umap.png"))
        markers = pd.read_csv(result_dir / "markers.csv")
        st.dataframe(markers)
        st.caption("Candidate annotations default to unknown. Export and review the table before reuse.")
        current = pd.read_csv(result_dir / "annotations.csv", keep_default_na=False)
        edited = st.data_editor(current, num_rows="fixed", hide_index=True)
        st.download_button("Download reviewed annotations CSV", edited.to_csv(index=False).encode("utf-8"), "annotations_reviewed.csv", "text/csv")
        st.caption("Pass the reviewed file to the CLI with --annotations for a new reproducible run.")
    else:
        st.info("Run an analysis to see results.")

with pages[3]:
    if (result_dir / "report.html").is_file():
        st.write(f"HTML report: {result_dir / 'report.html'}")
        for filename in ("report.html", "run_manifest.json", "qc_metrics.csv", "markers.csv", "annotations.csv", "processed.h5ad"):
            st.download_button(filename, (result_dir / filename).read_bytes(), filename, key=filename)
        st.caption("AI interpretation is disabled in v0.1; all analysis and reporting run locally without an API key.")
    else:
        st.info("No report in the selected results folder yet.")
