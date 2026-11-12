import hashlib
import json
from pathlib import Path
import zipfile

import pandas as pd
import pytest

from embryoomics.config import InputError
from embryoomics.evidence import load_annotations
from embryoomics.reporting import write_report
from embryoomics.review import annotation_version, save_annotation_review


@pytest.fixture
def review_run(tmp_path):
    run = tmp_path / 'review_run'
    run.mkdir()
    annotations = load_annotations(None, ['0', '1'])
    annotations.to_csv(run / 'annotations.csv', index=False)
    pd.DataFrame({'cluster': ['0', '1'], 'cells': [3, 4]}).to_csv(run / 'cluster_composition.csv', index=False)
    pd.DataFrame({'sample_id': ['test'], 'retained_cells': [7]}).to_csv(run / 'qc_preview.csv', index=False)
    pd.DataFrame({'cluster_id': ['0'], 'gene_symbol': ['Gene1']}).to_csv(run / 'markers.csv', index=False)
    (run / 'processed.h5ad').write_bytes(b'immutable analysis artifact')
    (run / 'figures').mkdir()
    for image in ('qc.png', 'umap.png'):
        (run / 'figures' / image).write_bytes(b'image placeholder')
    manifest = {'dataset_accession': 'SYNTHETIC', 'sample_id': 'test', 'stage': 'test',
                'source_url': 'https://example.invalid', 'input_dimensions': {'cells': 7, 'genes': 1},
                'genome_build': 'unknown', 'data_license': 'unknown', 'seed': 42, 'warnings': [],
                'files': ['report.html', 'run_manifest.json', 'annotations.csv', 'cluster_composition.csv',
                          'qc_preview.csv', 'markers.csv', 'processed.h5ad', 'figures/qc.png', 'figures/umap.png']}
    write_report(run, manifest, pd.read_csv(run / 'qc_preview.csv'), pd.read_csv(run / 'markers.csv'), annotations)
    return run, annotations


def snapshot(run):
    return {p.relative_to(run).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in run.rglob('*') if p.is_file()}


def test_annotation_save_refreshes_report_and_preserves_analysis(review_run):
    run, edited = review_run
    before = snapshot(run)
    version = annotation_version(run)
    edited.loc[0, 'candidate_cell_type'] = '候选 <未确认>'
    edited.loc[0, 'review_status'] = 'in_review'
    result = save_annotation_review(run, edited, version)
    assert result['changed']
    saved = pd.read_csv(run / 'annotations.csv', keep_default_na=False)
    assert saved.loc[0, 'candidate_cell_type'] == '候选 <未确认>'
    assert '候选 &lt;未确认&gt;' in (run / 'report.html').read_text(encoding='utf-8')
    for name in ('processed.h5ad', 'markers.csv', 'figures/umap.png', 'cluster_composition.csv'):
        assert hashlib.sha256((run / name).read_bytes()).hexdigest() == before[name]
    manifest = json.loads((run / 'run_manifest.json').read_text(encoding='utf-8'))
    history = json.loads((run / manifest['annotation_review']['history_file']).read_text(encoding='utf-8'))
    assert history['before'][0]['candidate_cell_type'] == 'unknown'
    assert history['after'][0]['candidate_cell_type'] == '候选 <未确认>'
    with zipfile.ZipFile(run / 'report_bundle.zip') as bundle:
        assert set(bundle.namelist()) == set(manifest['files'])
        assert bundle.read('annotations.csv') == (run / 'annotations.csv').read_bytes()
        assert bundle.testzip() is None
    with pytest.raises(InputError, match='CONFLICT'):
        save_annotation_review(run, edited, version)
    assert not save_annotation_review(run, edited, result['version'])['changed']


@pytest.mark.parametrize('column,value', [('cluster_id', '99'), ('candidate_cell_type', '  '), ('review_status', 'invalid')])
def test_invalid_edits_leave_run_unchanged(review_run, column, value):
    run, edited = review_run
    before = snapshot(run)
    edited.loc[0, column] = value
    with pytest.raises(InputError, match='ANNOTATIONS'):
        save_annotation_review(run, edited, annotation_version(run))
    assert snapshot(run) == before


def test_report_failure_does_not_partially_save(review_run, monkeypatch):
    import embryoomics.review as review
    run, edited = review_run
    before = snapshot(run)
    edited.loc[0, 'reviewer_note'] = 'draft'
    def fail(*args, **kwargs):
        raise OSError('report generation failed')
    monkeypatch.setattr(review, 'write_report', fail)
    with pytest.raises(OSError):
        save_annotation_review(run, edited, annotation_version(run))
    assert snapshot(run) == before


def test_commit_failure_rolls_back_saved_files(review_run, monkeypatch):
    import embryoomics.review as review
    run, edited = review_run
    before = snapshot(run)
    edited.loc[0, 'reviewer_note'] = 'draft'
    original_replace = review.os.replace
    def fail_manifest(source, target):
        if Path(target) == run / 'run_manifest.json' and Path(source).parent.name != 'backup':
            raise OSError('cannot replace manifest')
        return original_replace(source, target)
    monkeypatch.setattr(review.os, 'replace', fail_manifest)
    with pytest.raises(OSError):
        save_annotation_review(run, edited, annotation_version(run))
    assert snapshot(run) == before
