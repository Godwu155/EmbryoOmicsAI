import gzip
import json
import zipfile
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
    assert json.loads((out / 'run_status.json').read_text())['status'] == 'completed'
    assert manifest['marker_plots']['status'] == 'generated'
    assert 'figures/marker_dotplot.png' in manifest['files']
    assert 'figures/marker_heatmap.png' in manifest['files']
    panel = pd.read_csv(out / 'marker_plot_genes.csv')
    assert panel.gene_id.is_unique and len(panel) <= 30
    with zipfile.ZipFile(out / 'report_bundle.zip') as bundle:
        assert set(bundle.namelist()) == set(manifest['files'])
        assert bundle.testzip() is None
    report = (out / "report.html").read_text(encoding="utf-8")
    assert "id='content-en'" in report and "id='content-zh-CN'" in report
    assert "Language / 语言" in report and "setReportLanguage" in report
    assert "UMAP shows similarity" in report and "UMAP 显示本次分析中的相似性" in report
    assert "gene_id" in pd.read_csv(out / "markers.csv").columns
    assert "figures/marker_dotplot.png" in report and "Marker 平均表达热图" in report
    assert set(pd.read_csv(out / "annotations.csv").candidate_cell_type) == {"unknown"}
    processed = ad.read_h5ad(out / "processed.h5ad")
    assert sparse.issparse(processed.layers["counts"])
    np.testing.assert_array_equal(processed.layers["counts"].toarray(), values)
    with pytest.raises(InputError, match="OUTPUT_EXISTS"):
        run_analysis(config_path, folder, out)
    # Bind the workspace explicitly so this also works with pytest's default temp directory.
    with monkeypatch.context() as workspace:
        workspace.chdir(tmp_path)
        assert run_analysis(config_path, folder, out, overwrite=True)["cells"] == 30
    other_workspace = tmp_path / 'other_workspace'
    other_workspace.mkdir()
    with monkeypatch.context() as workspace:
        workspace.chdir(other_workspace)
        with pytest.raises(InputError, match='OVERWRITE_UNSAFE'):
            run_analysis(config_path, folder, out, overwrite=True)
    unsafe = tmp_path / "unrelated"
    unsafe.mkdir()
    with pytest.raises(InputError, match="OVERWRITE_UNSAFE"):
        run_analysis(config_path, folder, unsafe, overwrite=True)
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run(timeout=30)
    app.selectbox(key='language').set_value("简体中文").run(timeout=30)
    app.text_input(key='out_path').set_value(str(out)).run(timeout=30)
    assert not app.exception


def test_streamlit_language_switch():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run(timeout=30)
    assert [tab.label for tab in app.tabs] == ["Data", "QC", "Results", "Report"]
    app.selectbox(key='language').set_value("简体中文").run(timeout=30)
    assert [tab.label for tab in app.tabs] == ["数据", "质控", "结果", "报告"]
    assert not app.exception


def test_qc_confirmation_invalidated_when_config_changes(sample):
    from streamlit.testing.v1 import AppTest
    folder, config_path, config, _ = sample
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run()
    app.text_input(key='input_path').set_value(str(folder))
    app.text_input(key='config_path').set_value(str(config_path)).run()
    app.button(key='validate').click().run()
    app.checkbox(key='qc_confirm').check().run()
    assert app.checkbox(key='qc_confirm').value
    config['qc']['min_genes'] = 2
    config_path.write_text(yaml.safe_dump(config), encoding='utf-8')
    app.run()
    assert not app.exception
    assert not app.checkbox
    assert not app.session_state['qc_confirm']


def test_gene_by_cell_tsv_requires_layout_and_preserves_orientation(tmp_path):
    path = tmp_path / "counts.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("Gene\tcell_a\tcell_b\nENSG1\t2\t0\nENSG2\t1\t3\n")
    config = {"species": "homo_sapiens", "raw_count_source": "X", "sample_id": "test", "stage": "blastocyst"}
    with pytest.raises(InputError, match="TSV_LAYOUT"):
        load_counts(path, config)
    config["count_matrix_layout"] = "genes_by_cells_tsv"
    data, _ = load_counts(path, config)
    assert data.shape == (2, 2)
    assert data.obs_names.tolist() == ["cell_a", "cell_b"]
    assert data.var_names.tolist() == ["ENSG1", "ENSG2"]
    np.testing.assert_array_equal(data.layers["counts"].toarray(), [[2, 1], [0, 3]])
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("Gene\tcell_a\tcell_b\nENSG1\t2.5\t0\n")
    with pytest.raises(InputError, match="COUNTS_NOT_RAW"):
        load_counts(path, config)


def test_direct_annotation_save_in_interface(sample, tmp_path):
    from streamlit.testing.v1 import AppTest
    folder, config_path, _, _ = sample
    out = tmp_path / 'editable_run'
    run_analysis(config_path, folder, out)
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    app.text_input(key='out_path').set_value(str(out)).run(timeout=30)
    editor_key = next(key for key in app.session_state.filtered_state if key.startswith('annotation_editor:'))
    app.session_state[editor_key] = {'edited_rows': {0: {'candidate_cell_type': 'test candidate', 'review_status': 'in_review'}},
                                   'added_rows': [], 'deleted_rows': []}
    app.button(key='annotations_save').click().run(timeout=30)
    assert not app.exception
    saved = pd.read_csv(out / 'annotations.csv')
    assert saved.loc[0, 'candidate_cell_type'] == 'test candidate'
    assert 'test candidate' in (out / 'report.html').read_text(encoding='utf-8')
    assert any('Annotations saved' in element.value for element in app.success)
    # A fresh page must load the saved edits, rather than session-only state.
    fresh = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    fresh.text_input(key='out_path').set_value(str(out)).run(timeout=30)
    assert not fresh.exception
    assert fresh.dataframe[-1].value.iloc[0].candidate_cell_type == 'test candidate'
    from embryoomics.review import annotation_version, save_annotation_review
    changed = pd.read_csv(out / 'annotations.csv', dtype=str, keep_default_na=False)
    changed.loc[0, 'reviewer_note'] = 'saved from another page'
    save_annotation_review(out, changed, annotation_version(out))
    fresh.run(timeout=30)
    assert fresh.button(key='annotations_save').disabled
    fresh.button(key='annotations_reload').click().run(timeout=30)
    assert not fresh.exception
    assert not fresh.button(key='annotations_save').disabled
    assert fresh.dataframe[-1].value.iloc[0].reviewer_note == 'saved from another page'


def test_qc_interface_preview_matches_run_and_preserves_source_config(sample, tmp_path):
    from streamlit.testing.v1 import AppTest
    folder, config_path, _, values = sample
    original_config = config_path.read_bytes()
    out = tmp_path / 'ui_qc_run'
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    app.text_input(key='input_path').set_value(str(folder))
    app.text_input(key='config_path').set_value(str(config_path))
    app.text_input(key='out_path').set_value(str(out)).run(timeout=30)
    app.button(key='validate').click().run(timeout=30)
    app.checkbox(key='qc_confirm').check().run(timeout=30)
    app.number_input(key='qc_value_min_genes').set_value(20).run(timeout=30)
    assert not app.checkbox(key='qc_confirm').value
    assert app.button(key='run').disabled
    app.checkbox(key='qc_enable_max_pct_mt').check().run(timeout=30)
    assert app.checkbox(key='qc_confirm').disabled  # enabled filter still needs a number
    app.number_input(key='qc_value_max_pct_mt').set_value(20.0).run(timeout=30)
    expected = int((((values > 0).sum(axis=1) >= 20) & (values[:, 0] / values.sum(axis=1) <= .2)).sum())
    assert expected >= 4
    assert app.session_state['qc_current_preview']['retained'] == expected
    app.selectbox(key='language').set_value('简体中文').run(timeout=30)
    assert app.number_input(key='qc_value_min_genes').value == 20
    assert app.session_state['qc_current_preview']['retained'] == expected
    app.checkbox(key='qc_confirm').check().run(timeout=30)
    app.button(key='run').click().run(timeout=30)
    assert not app.exception
    assert config_path.read_bytes() == original_config
    manifest = json.loads((out / 'run_manifest.json').read_text(encoding='utf-8'),
                          parse_constant=lambda text: pytest.fail(f'Nonstandard JSON value {text}'))
    assert manifest['qc']['retained_cells'] == expected
    assert manifest['parameters']['qc'] == {'min_genes': 20, 'max_genes': None, 'max_pct_mt': 20.0}
    assert yaml.safe_load((out / 'config.yaml').read_text())['qc'] == manifest['parameters']['qc']
    assert pd.read_csv(out / 'qc_filter_reasons.csv').newly_removed_cells.sum() == 30 - expected
    processed = ad.read_h5ad(out / 'processed.h5ad')
    assert processed.n_obs == expected


def test_qc_interface_rejects_invalid_range_and_empty_result(sample):
    from streamlit.testing.v1 import AppTest
    folder, config_path, _, _ = sample
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    app.text_input(key='input_path').set_value(str(folder))
    app.text_input(key='config_path').set_value(str(config_path)).run(timeout=30)
    app.button(key='validate').click().run(timeout=30)
    app.number_input(key='qc_value_min_genes').set_value(99999).run(timeout=30)
    assert not app.exception
    assert app.button(key='run').disabled
    assert app.session_state['qc_current_preview']['retained'] == 0
    app.checkbox(key='qc_enable_max_genes').check().run(timeout=30)
    app.number_input(key='qc_value_max_genes').set_value(10).run(timeout=30)
    assert app.button(key='run').disabled
    assert any('min_genes exceeds' in item.value for item in app.error)
    app.checkbox(key='qc_enable_min_genes').uncheck().run(timeout=30)
    assert app.session_state['qc_current_preview']['thresholds']['min_genes'] is None


def test_stale_qc_preview_rejected_before_run(sample, tmp_path):
    folder, config_path, config, _ = sample
    info = inspect_input(folder, config_path)
    config['seed'] = 43
    config_path.write_text(yaml.safe_dump(config), encoding='utf-8')
    out = tmp_path / 'stale_run'
    with pytest.raises(InputError, match='QC_STALE'):
        run_analysis(config_path, folder, out, qc_override=info['config']['qc'], expected_input_identity=info['input_identity'])
    assert not out.exists()


def test_missing_mitochondrial_symbols_disable_ui_filter(sample):
    from streamlit.testing.v1 import AppTest
    folder, config_path, config, _ = sample
    feature_file = folder / 'GSMtest_features.tsv'
    feature_file.write_text(feature_file.read_text(encoding='utf-8').replace('mt-Nd1', 'Gene0'), encoding='utf-8')
    config['qc']['max_pct_mt'] = 10
    config_path.write_text(yaml.safe_dump(config), encoding='utf-8')
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    app.text_input(key='input_path').set_value(str(folder))
    app.text_input(key='config_path').set_value(str(config_path)).run(timeout=30)
    app.button(key='validate').click().run(timeout=30)
    assert not app.exception
    assert app.checkbox(key='qc_enable_max_pct_mt').disabled
    assert not app.checkbox(key='qc_enable_max_pct_mt').value
    assert app.session_state['qc_current_preview']['thresholds']['max_pct_mt'] is None
    assert app.session_state['qc_current_preview']['retained'] == 30


def test_result_explorer_and_selected_cluster_annotation_save(sample, tmp_path):
    from streamlit.testing.v1 import AppTest
    folder, config_path, _, _ = sample
    out = tmp_path / 'explorer_run'
    run_analysis(config_path, folder, out)
    before = pd.read_csv(out / 'annotations.csv', dtype=str, keep_default_na=False)
    cluster = before.cluster_id.iloc[0]
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    app.text_input(key='out_path').set_value(str(out)).run(timeout=30)
    app.selectbox(key=f'result_cluster:{out.resolve()}').set_value(cluster).run(timeout=30)
    markers = next(item.value for item in app.dataframe if 'gene_id' in item.value and 'cluster_id' in item.value)
    assert set(markers.cluster_id) == {cluster}
    assert app.dataframe[-1].value.cluster_id.tolist() == [cluster]
    app.text_input(key=f'gene_query:{out.resolve()}').set_value('Gene1').run(timeout=30)
    gene_key = next(key for key in app.session_state.filtered_state if key.startswith('gene_selection:'))
    app.selectbox(key=gene_key).set_value('id1').run(timeout=30)
    assert not app.exception
    summary = next(item.value for item in app.dataframe if 'mean_log1p_expression' in item.value)
    assert set(summary.cluster) == set(before.cluster_id)
    editor_key = next(key for key in app.session_state.filtered_state if key.startswith('annotation_editor:'))
    app.session_state[editor_key] = {'edited_rows': {0: {'reviewer_note': 'selected cluster note'}}, 'added_rows': [], 'deleted_rows': []}
    app.button(key='annotations_save').click().run(timeout=30)
    assert not app.exception
    after = pd.read_csv(out / 'annotations.csv', dtype=str, keep_default_na=False)
    assert after.loc[after.cluster_id == cluster, 'reviewer_note'].iloc[0] == 'selected cluster note'
    pd.testing.assert_frame_equal(before.loc[before.cluster_id != cluster].reset_index(drop=True),
                                  after.loc[after.cluster_id != cluster].reset_index(drop=True))
    app.text_input(key=f'gene_query:{out.resolve()}').set_value('not-a-gene').run(timeout=30)
    assert not app.exception
    assert any('No matching gene' in item.value for item in app.info)
