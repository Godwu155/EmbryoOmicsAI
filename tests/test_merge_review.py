import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
import yaml

from embryoomics.config import InputError
from embryoomics.io import load_counts
from embryoomics.pipeline import run_analysis
from embryoomics.review import annotation_snapshot, save_annotation_review
from test_annotation_review import review_run
from test_contract import sample


def test_annotation_snapshot_keeps_display_and_version_consistent(review_run, monkeypatch):
    run, _ = review_run
    target = run / 'annotations.csv'
    original = target.read_bytes()
    real_read = Path.read_bytes
    changed = False
    def concurrent_save(path):
        nonlocal changed
        content = real_read(path)
        if path == target and not changed:
            changed = True
            target.write_bytes(content.replace(b'unknown', b'new candidate'))
        return content
    monkeypatch.setattr(Path, 'read_bytes', concurrent_save)
    shown, version = annotation_snapshot(run)
    assert version == hashlib.sha256(original).hexdigest()
    assert set(shown.candidate_cell_type) == {'unknown'}
    shown.loc[0, 'reviewer_note'] = 'stale draft'
    with pytest.raises(InputError, match='ANNOTATIONS_CONFLICT'):
        save_annotation_review(run, shown, version)
    assert b'new candidate' in real_read(target)


@pytest.mark.parametrize('contained', ['counts', 'config'])
def test_overwrite_cannot_remove_source_files(sample, tmp_path, monkeypatch, contained):
    folder, config_path, config, values = sample
    out = tmp_path / 'old_run'
    out.mkdir()
    (out / 'run_manifest.json').write_text('{}')
    if contained == 'counts':
        data = ad.AnnData(sparse.csr_matrix(values))
        data.var['gene_symbol'] = [f'Gene{i}' for i in range(data.n_vars)]
        source = out / 'input.h5ad'
        data.write_h5ad(source)
    else:
        source = folder
        config_path = out / 'config.yaml'
        config_path.write_text(yaml.safe_dump(config))
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    monkeypatch.chdir(tmp_path)
    with pytest.raises(InputError, match='OVERWRITE_INPUT'):
        run_analysis(config_path, source, out, overwrite=True)
    assert {p.name: p.read_bytes() for p in out.iterdir()} == before


@pytest.mark.parametrize('dtype', [np.uint64, np.float64])
def test_h5ad_counts_outside_int64_rejected(sample, tmp_path, dtype):
    _, _, config, _ = sample
    data = ad.AnnData(sparse.csr_matrix(np.array([[2**63]], dtype=dtype)))
    path = tmp_path / 'overflow.h5ad'
    data.write_h5ad(path)
    with pytest.raises(InputError, match='COUNTS_RANGE'):
        load_counts(path, config)
