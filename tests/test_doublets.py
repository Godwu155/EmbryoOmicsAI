import json
from pathlib import Path
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
import yaml

from embryoomics.config import InputError, validate_config
from embryoomics.doublets import score_doublets, apply_doublet_decision
from embryoomics.pipeline import inspect_doublets, run_analysis
from embryoomics.review import annotation_version, save_annotation_review
from test_contract import sample


SETTINGS = {'mode': 'score', 'expected_doublet_rate': 0.05, 'threshold': 0.15,
            'n_prin_comps': 2, 'upstream_handling': 'unknown'}


@pytest.mark.parametrize('key,value', [('mode', 'auto'), ('expected_doublet_rate', None),
    ('threshold', float('nan')), ('threshold', True), ('expected_doublet_rate', 1),
    ('n_prin_comps', True), ('n_prin_comps', 0), ('upstream_handling', 'guessed')])
def test_invalid_doublet_config(sample, key, value):
    config = sample[2]
    config['doublets'] = {**SETTINGS, key: value}
    with pytest.raises(InputError, match='CONFIG'):
        validate_config(config)


def test_scrublet_raw_counts_and_reproducibility():
    rng = np.random.default_rng(34)
    counts = rng.poisson(1, (120, 180))
    counts[:60, :30] += rng.poisson(4, (60, 30))
    counts[60:, 30:60] += rng.poisson(4, (60, 30))
    data = ad.AnnData(sparse.csr_matrix(np.ones_like(counts)))
    data.obs['sample_id'] = 'sample'
    data.layers['counts'] = sparse.csr_matrix(counts)
    scores, simulated, summary = score_doublets(data, SETTINGS, 42)
    raw = ad.AnnData(sparse.csr_matrix(counts), obs=data.obs.copy())
    raw.layers['counts'] = raw.X.copy()
    repeat, _, repeated = score_doublets(raw, SETTINGS, 42)
    pd.testing.assert_frame_equal(scores, repeat)
    assert summary['score_digest'] == repeated['score_digest']
    assert len(scores) == 120 and np.isfinite(simulated).all()
    np.testing.assert_array_equal(data.X.toarray(), np.ones_like(counts))
    np.testing.assert_array_equal(data.layers['counts'].toarray(), counts)
    data.obs['batch'] = ['a'] * 60 + ['b'] * 60
    with pytest.raises(InputError, match='SCRUBLET_BATCH'):
        score_doublets(data, SETTINGS, 42)
    with pytest.raises(InputError, match='SCRUBLET_SIZE'):
        score_doublets(raw[:5], SETTINGS, 42)


def test_qc_failed_cells_not_scored_and_filtering_is_explicit():
    metrics = pd.DataFrame(index=['a', 'b', 'c'])
    scores = pd.DataFrame({'doublet_score': [.1, .3], 'predicted_doublet': [False, True]}, index=['a', 'c'])
    mask, audit = apply_doublet_decision(metrics, [True, False, True], scores, SETTINGS)
    assert mask.tolist() == [True, False, True]
    assert pd.isna(audit.loc['b', 'doublet_score']) and not audit.removed_as_doublet.any()
    mask, audit = apply_doublet_decision(metrics, [True, False, True], scores, {**SETTINGS, 'mode': 'filter'})
    assert mask.tolist() == [True, False, False]
    assert audit.loc['c', 'removed_as_doublet']


def test_pipeline_preview_filter_and_report_audit(sample, tmp_path):
    folder, config_path, config, values = sample
    config['doublets'] = {**SETTINGS, 'mode': 'filter'}
    config_path.write_text(yaml.safe_dump(config), encoding='utf-8')
    out = tmp_path / 'run'
    with pytest.raises(InputError, match='DOUBLET_CONFIRMATION'):
        run_analysis(config_path, folder, out)
    preview = inspect_doublets(folder, config_path)
    with pytest.raises(InputError, match='DOUBLET_STALE'):
        run_analysis(config_path, folder, out, confirm_doublets=True, expected_doublet_digest='old')
    run_analysis(config_path, folder, out, confirm_doublets=True,
                 expected_input_identity=preview['input_identity'], expected_doublet_digest=preview['summary']['score_digest'])
    manifest = json.loads((out / 'run_manifest.json').read_text(encoding='utf-8'))
    audit = pd.read_csv(out / 'doublet_scores.csv', index_col='cell_id')
    assert manifest['doublets']['filtering_confirmed']
    assert int(audit.retained.sum()) == manifest['qc']['retained_cells'] == preview['retained']
    processed = ad.read_h5ad(out / 'processed.h5ad')
    np.testing.assert_array_equal(processed.layers['counts'].toarray(), values[audit.retained.to_numpy()])
    assert processed.obs_names.tolist() == audit.index[audit.retained].tolist()
    assert not processed.obs.predicted_doublet.any()
    report = (out / 'report.html').read_text(encoding='utf-8')
    assert 'Scrublet 审核与过滤决定' in report and 'figures/doublet_scores.png' in report
    with zipfile.ZipFile(out / 'report_bundle.zip') as bundle:
        assert set(bundle.namelist()) == set(manifest['files']) and bundle.testzip() is None
    annotations = pd.read_csv(out / 'annotations.csv', dtype=str, keep_default_na=False)
    annotations.loc[0, 'reviewer_note'] = 'software audit'
    save_annotation_review(out, annotations, annotation_version(out))
    assert 'figures/doublet_scores.png' in (out / 'report.html').read_text(encoding='utf-8')


def test_score_only_keeps_candidates(sample, tmp_path):
    folder, config_path, config, values = sample
    config['doublets'] = SETTINGS
    config['qc']['min_genes'] = 18
    config_path.write_text(yaml.safe_dump(config), encoding='utf-8')
    result = run_analysis(config_path, folder, tmp_path / 'score')
    expected = int(((values > 0).sum(axis=1) >= 18).sum())
    assert 10 <= expected < 30
    assert result['cells'] == expected
    audit = pd.read_csv(tmp_path / 'score/doublet_scores.csv')
    assert audit.doublet_score.isna().sum() == 30 - expected
    assert ad.read_h5ad(tmp_path / 'score/processed.h5ad').obs.predicted_doublet.dtype == bool
    summary = json.loads((tmp_path / 'score/run_manifest.json').read_text(encoding='utf-8'))['doublets']
    assert summary['removed_cells'] == 0 and not summary['filtering_confirmed']


def test_cli_preview_exports_and_requires_separate_filter_confirmation(sample, tmp_path, capsys):
    from embryoomics.cli import main
    folder, config_path, config, _ = sample
    config['doublets'] = {**SETTINGS, 'mode': 'filter'}
    config_path.write_text(yaml.safe_dump(config), encoding='utf-8')
    out = tmp_path / 'cli-preview'
    assert main(['preview-doublets', '--input', str(folder), '--config', str(config_path), '--out', str(out)]) == 0
    saved = json.loads((out / 'doublet_preview.json').read_text(encoding='utf-8'))
    assert saved['config']['doublets']['mode'] == 'filter'
    assert len(pd.read_csv(out / 'doublet_scores.csv')) == 30
    assert main(['run', '--input', str(folder), '--config', str(config_path), '--out', str(tmp_path / 'run'), '--confirm-qc']) == 2
    assert 'DOUBLET_CONFIRMATION' in capsys.readouterr().err


def test_interface_requires_preview_confirmation_and_invalidates_changes(sample, tmp_path):
    from streamlit.testing.v1 import AppTest
    folder, config_path, _, _ = sample
    source = config_path.read_bytes()
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'app.py')).run(timeout=30)
    app.text_input(key='input_path').set_value(str(folder))
    app.text_input(key='config_path').set_value(str(config_path))
    out = tmp_path / 'ui_run'
    app.text_input(key='out_path').set_value(str(out)).run(timeout=30)
    app.button(key='validate').click().run(timeout=30)
    app.selectbox(key='doublet_mode').set_value('filter').run(timeout=30)
    assert app.button(key='run').disabled
    app.number_input(key='doublet_rate').set_value(0.05)
    app.number_input(key='doublet_threshold').set_value(0.15)
    app.number_input(key='doublet_pcs').set_value(2).run(timeout=30)
    assert app.button(key='run').disabled
    app.button(key='doublet_preview_button').click().run(timeout=30)
    assert not app.exception
    assert app.checkbox(key='qc_confirm').disabled
    app.checkbox(key='doublet_confirm').check().run(timeout=30)
    app.checkbox(key='qc_confirm').check().run(timeout=30)
    assert not app.button(key='run').disabled
    app.number_input(key='doublet_threshold').set_value(0.16).run(timeout=30)
    assert app.button(key='run').disabled
    assert not app.checkbox(key='doublet_confirm').value
    assert 'doublet_preview' not in app.session_state
    app.button(key='doublet_preview_button').click().run(timeout=30)
    app.checkbox(key='doublet_confirm').check().run(timeout=30)
    app.checkbox(key='qc_confirm').check().run(timeout=30)
    expected = app.session_state['doublet_preview']['retained']
    app.selectbox(key='language').set_value('简体中文').run(timeout=30)
    assert app.number_input(key='doublet_threshold').value == 0.16
    app.button(key='run').click().run(timeout=30)
    assert not app.exception
    assert config_path.read_bytes() == source
    manifest = json.loads((out / 'run_manifest.json').read_text(encoding='utf-8'))
    assert manifest['qc']['retained_cells'] == expected
    assert manifest['doublets']['filtering_confirmed']
