import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from embryoomics.cellxgene import export_view, import_proposals, apply_proposals, sha256
from embryoomics.config import InputError
from embryoomics.evidence import load_annotations
from test_annotation_review import review_run


@pytest.fixture
def exchange(tmp_path):
    run = tmp_path / 'run'
    run.mkdir()
    data = ad.AnnData(sparse.csr_matrix([[0., 1., 2.], [3., 0., 1.], [1., 2., 0.], [2., 1., 3.]]))
    data.obs_names = ['a', 'b', 'c', 'd']
    data.var_names = ['id1', 'id2', 'id3']
    data.var['gene_symbol'] = pd.Categorical(['Same', 'Same', None])
    data.obs['cluster'] = pd.Categorical(['0', '0', '1', '1'])
    data.obs['stage'], data.obs['sample_id'] = 'test', 'test'
    data.layers['counts'] = sparse.csr_matrix(np.ones((4, 3), dtype=int) * 10)
    data.obsm['X_umap'] = np.arange(8).reshape(4, 2).astype(float)
    data.write_h5ad(run / 'processed.h5ad')
    annotations = load_annotations(None, ['0', '1'])
    annotations.loc[0, 'support_genes'] = 'existing evidence'
    annotations.to_csv(run / 'annotations.csv', index=False)
    out = tmp_path / 'exchange'
    export_view(run, out)
    return run, out, data


def test_export_preserves_source_expression_embedding_and_distinct_ids(exchange):
    run, out, original = exchange
    view = ad.read_h5ad(out / 'viewer.h5ad')
    assert sparse.isspmatrix_csc(view.X) and view.X.dtype == np.float32
    np.testing.assert_allclose(view.X.toarray(), original.X.toarray())
    np.testing.assert_array_equal(view.obsm['X_umap'], original.obsm['X_umap'])
    assert view.var.viewer_gene.tolist() == ['Same [id1]', 'Same [id2]', 'id3 [id3]']
    assert not view.layers and view.raw is None
    provenance = json.loads((out / 'exchange.json').read_text(encoding='utf-8'))
    assert provenance['processed_sha256'] == sha256(run / 'processed.h5ad')
    assert provenance['annotations_sha256'] == sha256(run / 'annotations.csv')
    with pytest.raises(InputError, match='OUTPUT_EXISTS'):
        export_view(run, out)


def test_unanimous_and_mixed_labels_propose_without_overwriting(exchange, tmp_path):
    run, out, _ = exchange
    before = (run / 'annotations.csv').read_bytes()
    labels = out / 'cell_annotations.csv'
    labels.write_text('# CELLxGENE generated header\ncell_id,candidate_cell_type\nd,B\nb,Type#1\na,Type#1\nc,A\n', encoding='utf-8')
    dest = tmp_path / 'proposals'
    result = import_proposals(out, labels, dest)
    proposals = pd.read_csv(dest / 'annotation_proposals.csv', dtype=str, keep_default_na=False).set_index('cluster_id')
    assert result['unresolved_clusters'] == 1
    assert proposals.loc['0', 'candidate_cell_type'] == 'Type#1'
    assert proposals.loc['0', 'support_genes'] == 'existing evidence'
    assert proposals.loc['1', 'candidate_cell_type'] == 'unknown'
    assert set(proposals.review_status) == {'in_review'}
    assert (run / 'annotations.csv').read_bytes() == before


@pytest.mark.parametrize('rows', ['a,A\na,A\nc,A\nd,A\n', 'a,A\nb,A\nc,A\n', 'a,A\nb,A\nc,A\nx,A\n'])
def test_proposals_reject_duplicate_missing_or_extra_cells(exchange, tmp_path, rows):
    _, out, _ = exchange
    labels = out / 'cell_annotations.csv'
    labels.write_text('cell_id,candidate_cell_type\n' + rows, encoding='utf-8')
    with pytest.raises(InputError, match='VIEWER_LABELS'):
        import_proposals(out, labels, tmp_path / 'invalid')


def test_source_annotation_change_blocks_stale_import(exchange, tmp_path):
    run, out, _ = exchange
    with (run / 'annotations.csv').open('a', encoding='utf-8') as handle:
        handle.write('\n')
    with pytest.raises(InputError, match='VIEWER_STALE'):
        import_proposals(out, out / 'cell_annotations.csv', tmp_path / 'stale')


def test_changed_viewer_is_rejected(exchange, tmp_path):
    _, out, _ = exchange
    with (out / 'viewer.h5ad').open('ab') as handle:
        handle.write(b'changed')
    with pytest.raises(InputError, match='VIEWER_STALE'):
        import_proposals(out, out / 'cell_annotations.csv', tmp_path / 'stale')


def test_reviewed_proposals_use_existing_history_and_stale_guard(review_run, tmp_path):
    run, annotations = review_run
    proposals = tmp_path / 'proposals'
    proposals.mkdir()
    manifest = {'source_run': str(run), 'source_processed_sha256': sha256(run / 'processed.h5ad'),
                'source_annotations_sha256': sha256(run / 'annotations.csv')}
    (proposals / 'proposal_manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    annotations.loc[0, 'candidate_cell_type'] = 'reviewed test candidate'
    annotations.loc[0, 'review_status'] = 'in_review'
    annotations.to_csv(proposals / 'annotation_proposals.csv', index=False)
    with pytest.raises(InputError, match='VIEWER_CONFIRMATION'):
        apply_proposals(proposals)
    result = apply_proposals(proposals, confirm_reviewed=True)
    assert result['changed'] and list((run / 'annotation_history').glob('*.json'))
    assert 'reviewed test candidate' in (run / 'report.html').read_text(encoding='utf-8')
    with pytest.raises(InputError, match='ANNOTATIONS_CONFLICT'):
        apply_proposals(proposals, confirm_reviewed=True)
