"""Local Streamlit review UI using the same analysis service as CLI."""

from pathlib import Path

import pandas as pd
import numpy as np
import streamlit as st
import yaml
import json

from embryoomics.config import InputError, load_config, validate_config
from embryoomics.i18n import LANGUAGES, columns_for, t, warning_for
from embryoomics.pipeline import inspect_input, inspect_doublets, run_analysis
from embryoomics.io import review_identity
from embryoomics.review import annotation_snapshot, save_annotation_review, REVIEW_STATUSES
from embryoomics.qc import preview_qc_metrics
from embryoomics.doublets import apply_doublet_decision
from embryoomics.explorer import result_metadata, find_genes, gene_expression, expression_summary


@st.cache_data(show_spinner=False, max_entries=3)
def cached_result_metadata(path, fingerprint):
    return result_metadata(path)


@st.cache_data(show_spinner=False, max_entries=8)
def cached_gene_expression(path, fingerprint, gene_id):
    return gene_expression(path, gene_id)


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

try:
    current_identity = review_identity(input_path, config, metadata_path or None)
except (InputError, OSError):
    current_identity = None
if st.session_state.get('review_identity') != current_identity:
    st.session_state.pop('inspection', None)
    st.session_state['qc_confirm'] = False

with pages[0]:
    fields = ("dataset_accession", "sample_id", "stage", "species", "genome_build", "source_url", "data_license", "raw_count_source")
    labels = columns_for(language)
    st.dataframe(pd.DataFrame([{t(language, "field"): labels.get(key, key),
                                t(language, "value"): t(language, "unknown") if config[key] == "unknown" else config[key]} for key in fields]), hide_index=True)
    if st.button(t(language, "validate"), key="validate"):
        try:
            st.session_state["inspection"] = inspect_input(input_path, config_path, metadata_path or None, apply_config_filters=False)
            st.session_state["inspected_paths"] = (input_path, config_path, metadata_path)
            st.session_state['review_identity'] = current_identity
            st.session_state['qc_confirm'] = False
        except (InputError, FileNotFoundError, ValueError) as exc:
            st.error(str(exc))
    if "inspection" in st.session_state:
        info = st.session_state["inspection"]
        st.success(f"{info['cells']} {t(language, 'cells')} × {info['genes']} {t(language, 'genes')}")
        mt_method = t(language, 'mt_configured' if info['mt']['method'] == 'configured mitochondrial gene IDs' else "mt_symbol_method" if info["mt"]["matched_genes"] else "mt_no_match")
        st.write(t(language, "mt_detection"), {t(language, "mt_method"): mt_method,
                                               t(language, "mt_matched"): info["mt"]["matched_genes"]})
        for warning in info["warnings"]:
            st.warning(warning_for(language, warning))

with pages[1]:
    st.write(t(language, "qc_intro"))
    if "inspection" in st.session_state and st.session_state.get("inspected_paths") == (input_path, config_path, metadata_path):
        info = st.session_state["inspection"]
        metrics = info["qc_metrics"]
        if st.session_state.get('qc_controls_identity') != current_identity:
            st.session_state.pop('doublet_controls', None)
            for key in ('doublet_mode', 'doublet_rate', 'doublet_threshold', 'doublet_pcs', 'doublet_upstream'):
                st.session_state.pop(key, None)
            for rule in ('min_genes', 'max_genes', 'max_pct_mt'):
                st.session_state.pop(f'qc_enable_{rule}', None)
                st.session_state.pop(f'qc_value_{rule}', None)
            st.session_state['qc_controls_identity'] = current_identity
        thresholds = {}
        mt_available = np.isfinite(metrics.pct_counts_mt.to_numpy()).any()
        for rule in ('min_genes', 'max_genes', 'max_pct_mt'):
            available = rule != 'max_pct_mt' or mt_available
            enabled = st.checkbox(f"{t(language, 'qc_enable')} {t(language, 'qc_' + rule)}",
                                  value=config['qc'][rule] is not None and available,
                                  disabled=not available, key=f'qc_enable_{rule}')
            is_percentage = rule == 'max_pct_mt'
            initial_value = config['qc'][rule]
            if is_percentage and initial_value is not None:
                initial_value = float(initial_value)
            value = st.number_input(t(language, 'qc_' + rule),
                                    min_value=0.0 if is_percentage else 0,
                                    max_value=100.0 if is_percentage else None,
                                    value=initial_value, step=0.5 if is_percentage else 1,
                                    disabled=not enabled or not available,
                                    placeholder=t(language, 'qc_enter_value'), key=f'qc_value_{rule}')
            thresholds[rule] = value if enabled and available else None
        if not mt_available:
            st.info(t(language, 'qc_mt_unavailable'))
        st.subheader(t(language, 'doublets_title'))
        st.caption(t(language, 'doublets_help'))
        doublets = {**config['doublets'], **st.session_state.get('doublet_controls', {})}
        doublets['mode'] = st.selectbox(t(language, 'doublets_mode'), ['off', 'score', 'filter'],
            index=['off', 'score', 'filter'].index(doublets['mode']), format_func=lambda v: t(language, 'doublets_' + v), key='doublet_mode')
        if doublets['mode'] != 'off':
            doublets['expected_doublet_rate'] = st.number_input(t(language, 'doublets_rate'), min_value=0.001, max_value=0.999,
                value=float(doublets['expected_doublet_rate']) if doublets['expected_doublet_rate'] is not None else None, step=0.01, key='doublet_rate')
            doublets['threshold'] = st.number_input(t(language, 'doublets_threshold'), min_value=0.001, max_value=0.999,
                value=float(doublets['threshold']) if doublets['threshold'] is not None else None, step=0.01, key='doublet_threshold')
            doublets['n_prin_comps'] = st.number_input(t(language, 'doublets_pcs'), min_value=2, value=doublets['n_prin_comps'], step=1, key='doublet_pcs')
            doublets['upstream_handling'] = st.selectbox(t(language, 'doublets_upstream'), ['unknown', 'not_removed', 'removed'],
                index=['unknown', 'not_removed', 'removed'].index(doublets['upstream_handling']),
                format_func=lambda v: t(language, 'doublets_upstream_' + v), key='doublet_upstream')
        st.session_state['doublet_controls'] = doublets.copy()
        preview_identity = (current_identity, tuple(thresholds.items()), tuple(doublets.items()))
        if preview_identity != st.session_state.get('qc_preview_identity'):
            st.session_state['qc_confirm'] = False
            st.session_state['doublet_confirm'] = False
            st.session_state.pop('doublet_preview', None)
        st.session_state['qc_preview_identity'] = preview_identity
        can_run = False
        try:
            for rule in thresholds:
                if st.session_state.get(f'qc_enable_{rule}') and thresholds[rule] is None:
                    raise InputError(t(language, 'qc_value_required'))
            effective_config = validate_config({**config, 'qc': thresholds, 'doublets': doublets})
            mask, summary, reasons = preview_qc_metrics(metrics, effective_config['qc'])
            doublet_ready = True
            if doublets['mode'] != 'off':
                if st.button(t(language, 'doublets_preview'), key='doublet_preview_button'):
                    st.session_state['doublet_preview'] = inspect_doublets(input_path, config_path, metadata_path or None,
                        qc_override=thresholds, doublet_override=doublets)
                preview = st.session_state.get('doublet_preview')
                doublet_ready = preview is not None
                if preview is None:
                    st.info(t(language, 'doublets_preview_required'))
                else:
                    ordinary_retained = int(mask.sum())
                    st.write({t(language, 'doublets_candidates'): preview['summary']['candidate_doublets'],
                              t(language, 'qc_retained_cells'): preview['retained']})
                    st.dataframe(preview['audit'].rename(columns=labels), hide_index=False)
                    for key, values in [('doublets_observed', preview['scores'].doublet_score), ('doublets_simulated', preview['simulated'])]:
                        counts, edges = np.histogram(values, bins=30)
                        st.caption(t(language, key))
                        st.bar_chart(pd.DataFrame({t(language, 'chart_count'): counts}, index=(edges[:-1] + edges[1:]) / 2))
                    st.download_button(t(language, 'doublets_download'), preview['audit'].to_csv(), 'doublet_scores.csv', 'text/csv')
                    if doublets['mode'] == 'filter':
                        mask, _ = apply_doublet_decision(metrics, mask, preview['scores'], doublets)
                        summary['retained_cells'] = int(mask.sum())
                        summary['removed_cells'] = summary.input_cells - summary.retained_cells
                    reasons = pd.DataFrame([*reasons.to_dict(orient='records'), {
                        'filter_rule': 'doublets', 'threshold': doublets['threshold'], 'enabled': doublets['mode'] == 'filter',
                        'failed_cells': preview['summary']['candidate_doublets'], 'newly_removed_cells': ordinary_retained - int(mask.sum()),
                        'retained_after_rule': int(mask.sum())}])
            st.session_state['qc_current_preview'] = {'thresholds': thresholds.copy(), 'retained': int(mask.sum())}
            cards = st.columns(3)
            for card, label, count in zip(cards, ('input_cells', 'retained_cells', 'removed_cells'),
                                          (len(metrics), int(mask.sum()), int((~mask).sum()))):
                card.metric(t(language, 'qc_' + label), count)
            st.dataframe(summary.rename(columns=labels), hide_index=True)
            shown_reasons = reasons.copy()
            shown_reasons['filter_rule'] = shown_reasons.filter_rule.map(lambda rule: t(language, 'qc_' + rule))
            st.dataframe(shown_reasons.rename(columns=labels), hide_index=True)
            st.caption(t(language, 'qc_reasons_help'))
            st.download_button(t(language, 'qc_download_config'), yaml.safe_dump(effective_config, sort_keys=False, allow_unicode=True),
                               'run_config.yaml', 'text/yaml', key='qc_config_download')
            if mask.sum() < 4:
                st.error(t(language, 'qc_too_few'))
            elif (metrics.loc[mask, 'total_counts'] == 0).any():
                st.error(t(language, 'qc_empty_cells'))
            else:
                can_run = doublet_ready
        except InputError as exc:
            st.session_state.pop('qc_current_preview', None)
            st.error(str(exc))
        if not can_run:
            st.session_state['qc_confirm'] = False
        for key in ("total_counts", "n_genes_by_counts", "pct_counts_mt"):
            if metrics[key].notna().any():
                st.subheader(labels.get(key, key))
                counts, edges = np.histogram(metrics[key].dropna().to_numpy(dtype=float), bins=30)
                histogram = pd.DataFrame({"bin_midpoint": (edges[:-1] + edges[1:]) / 2, t(language, "chart_count"): counts}).set_index("bin_midpoint")
                st.bar_chart(histogram)
        st.caption(t(language, "qc_edit"))
        doublet_confirmed = False
        if doublets['mode'] == 'filter':
            doublet_confirmed = st.checkbox(t(language, 'doublets_confirm'), disabled=not can_run, key='doublet_confirm')
            can_run = can_run and doublet_confirmed
        confirmed = st.checkbox(t(language, "qc_confirm"), disabled=not can_run, key="qc_confirm")
        if st.button(t(language, "run"), disabled=not confirmed or not can_run, key="run"):
            try:
                result = run_analysis(config_path, input_path, out_path, metadata_path or None,
                                      qc_override=thresholds, expected_input_identity=info['input_identity'], doublet_override=doublets,
                                      confirm_doublets=doublet_confirmed,
                                      expected_doublet_digest=st.session_state['doublet_preview']['summary']['score_digest'] if doublets['mode'] != 'off' else None)
                st.session_state["result"] = result
                st.success(f"{t(language, 'completed')}: {result['output']}")
            except (InputError, FileNotFoundError, ValueError) as exc:
                st.error(str(exc))
    else:
        st.info(t(language, "validate_first"))

result_dir = Path(out_path)
with pages[2]:
    if (result_dir / 'run_status.json').is_file():
        run_status = json.loads((result_dir / 'run_status.json').read_text(encoding='utf-8'))
        if run_status.get('status') != 'completed':
            st.error(t(language, 'run_incomplete'))
            st.code(run_status.get('message', ''))
            st.stop()
    if (result_dir / 'figures/doublet_scores.png').is_file():
        with st.expander(t(language, 'doublets_summary')):
            st.image(str(result_dir / 'figures/doublet_scores.png'))
            st.caption(t(language, 'doublets_help'))
            st.download_button(t(language, 'doublets_download'), (result_dir / 'doublet_scores.csv').read_bytes(),
                               'doublet_scores.csv', 'text/csv', key='saved_doublet_download')
    if (result_dir / "markers.csv").is_file():
        current, version = annotation_snapshot(result_dir)
        selected_cluster = st.selectbox(t(language, 'result_cluster'), ['__all__', *current.cluster_id.tolist()],
                                        format_func=lambda value: t(language, 'all_clusters') if value == '__all__' else value,
                                        key=f'result_cluster:{result_dir.resolve()}')
        st.caption(f"{t(language, 'umap_axes')} {t(language, 'umap_limit')}")
        try:
            processed = result_dir / 'processed.h5ad'
            stat = processed.stat()
            fingerprint = (stat.st_size, stat.st_mtime_ns)
            points, genes = cached_result_metadata(str(processed.resolve()), fingerprint)
            selected = np.ones(len(points), dtype=bool) if selected_cluster == '__all__' else points.cluster.to_numpy() == selected_cluster
            st.scatter_chart(points.loc[selected], x='UMAP1', y='UMAP2', color='cluster', size=12, height=350)
            query = st.text_input(t(language, 'gene_search'), key=f'gene_query:{result_dir.resolve()}')
            if query.strip():
                matches = find_genes(genes, query)
                if matches.empty:
                    st.info(t(language, 'gene_no_match'))
                else:
                    st.caption(t(language, 'gene_search_limit'))
                    gene_id = st.selectbox(t(language, 'gene_select'), matches.index.tolist(), index=None,
                                           format_func=lambda gene: f"{genes.loc[gene, 'gene_symbol'] or gene} ({gene})",
                                           key=f'gene_selection:{result_dir.resolve()}:{query.strip().lower()}')
                    if gene_id is not None:
                        values = cached_gene_expression(str(processed.resolve()), fingerprint, gene_id)
                        gene_points = points.assign(log1p_expression=values)
                        st.scatter_chart(gene_points.loc[selected], x='UMAP1', y='UMAP2', color='log1p_expression', size=12, height=350)
                        st.caption(t(language, 'gene_expression_help'))
                        st.dataframe(expression_summary(points, values).rename(columns=labels), hide_index=True)
        except (InputError, OSError, KeyError, ValueError) as exc:
            st.warning(f"{t(language, 'result_view_error')}: {exc}")
            if (result_dir / 'figures' / 'umap.png').is_file():
                st.image(str(result_dir / 'figures' / 'umap.png'))
        for name in ('marker_dotplot', 'marker_heatmap'):
            figure = result_dir / 'figures' / f'{name}.png'
            if figure.is_file():
                with st.expander(t(language, name)):
                    st.image(str(figure))
                    st.caption(t(language, 'marker_plots_help'))
        markers = pd.read_csv(result_dir / "markers.csv", dtype={'cluster_id': str})
        if selected_cluster != '__all__':
            markers = markers.loc[markers.cluster_id == selected_cluster]
        shown_markers = markers.copy()
        if language == "zh-CN" and "test_method" in shown_markers:
            shown_markers["test_method"] = shown_markers["test_method"].replace({
                "Scanpy Wilcoxon on log1p normalized expression; cluster vs rest": "Scanpy Wilcoxon 检验：log1p 归一化表达，聚类与其余细胞比较"})
        st.dataframe(shown_markers.rename(columns=labels))
        st.caption(t(language, "annotation_hint"))
        st.caption(f"{t(language, 'annotation_target')}: {result_dir.resolve()}")
        shown_annotations = current if selected_cluster == '__all__' else current.loc[current.cluster_id == selected_cluster]
        st.caption(t(language, 'annotation_switch_hint'))
        context_key = f'annotation_version:{result_dir.resolve()}'
        expected_version = st.session_state.setdefault(context_key, version)
        stale = expected_version != version
        if stale:
            st.warning(t(language, 'annotation_conflict'))
            if st.button(t(language, 'annotation_reload'), key='annotations_reload'):
                st.session_state[context_key] = version
                st.rerun()
        flash = st.session_state.pop('annotation_saved', None)
        if flash:
            st.success(t(language, 'annotation_saved' if flash == 'changed' else 'annotation_unchanged'))
        with st.form('annotation_review_form'):
            columns = {key: st.column_config.TextColumn(label) for key, label in labels.items() if key in current}
            columns['review_status'] = st.column_config.SelectboxColumn(labels.get('review_status', 'Review status'), options=REVIEW_STATUSES, required=True)
            edited = st.data_editor(shown_annotations, num_rows="fixed", hide_index=True, disabled=['cluster_id'],
                                    column_config=columns, key=f'annotation_editor:{result_dir.resolve()}:{expected_version}:{selected_cluster}')
            save_clicked = st.form_submit_button(t(language, 'annotation_save'), disabled=stale, key='annotations_save')
        merged = pd.concat([current.loc[~current.cluster_id.isin(edited.cluster_id)], edited], ignore_index=True)
        if save_clicked:
            try:
                saved = save_annotation_review(result_dir, merged, expected_version)
                st.session_state[context_key] = saved['version']
                st.session_state['annotation_saved'] = 'changed' if saved['changed'] else 'unchanged'
                st.rerun()
            except (InputError, OSError, ValueError) as exc:
                st.error(str(exc))
        st.download_button(t(language, "annotation_download"), merged.to_csv(index=False).encode("utf-8"), "annotations_reviewed.csv", "text/csv", key="annotations_download")
        st.caption(t(language, "annotation_rerun"))
    else:
        st.info(t(language, "run_first"))

with pages[3]:
    if (result_dir / "report.html").is_file():
        st.write(f"{t(language, 'report_location')}: {result_dir / 'report.html'}")
        if (result_dir / 'report_bundle.zip').is_file():
            st.download_button('完整离线报告 / Complete offline report (ZIP)',
                               (result_dir / 'report_bundle.zip').read_bytes(), 'report_bundle.zip',
                               'application/zip', key='report_bundle')
        for filename in ("report.html", "run_manifest.json", "qc_metrics.csv", "markers.csv", "annotations.csv", "processed.h5ad"):
            st.download_button(filename, (result_dir / filename).read_bytes(), filename, key=filename)
        st.caption(t(language, "ai_disabled"))
    else:
        st.info(t(language, "report_first"))
