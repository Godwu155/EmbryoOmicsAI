import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from embryoomics.config import InputError
from embryoomics.explorer import result_metadata, find_genes, gene_expression, expression_summary


@pytest.mark.parametrize('format', ['csr', 'csc', 'dense'])
def test_read_selected_gene_and_metadata(tmp_path, format, monkeypatch):
    import embryoomics.explorer as explorer
    values = np.array([[0., 2., 0.], [3., 0., 1.], [0., 0., 0.], [4., 5., 6.]])
    matrix = sparse.csr_matrix(values) if format == 'csr' else sparse.csc_matrix(values) if format == 'csc' else values
    data = ad.AnnData(matrix)
    data.obs_names = ['c0', 'c1', 'c2', 'c3']
    data.var_names = ['id0', 'id1', 'id2']
    data.obs['cluster'] = ['0', '0', '1', '1']
    data.obs['sample_id'], data.obs['stage'] = 'test', 'test'
    data.var['gene_symbol'] = ['GeneA', 'GeneA', 'GeneB']
    data.obsm['X_umap'] = np.arange(8).reshape(4, 2).astype(float)
    path = tmp_path / 'processed.h5ad'
    data.write_h5ad(path)
    # Exercise chunk boundaries and empty sparse rows without using a full-matrix reader.
    monkeypatch.setattr(explorer, '_SPARSE_CHUNK_SIZE', 2)
    def forbid_full_read(*args, **kwargs):
        pytest.fail('Explorer must not load the complete AnnData object')
    monkeypatch.setattr(ad, 'read_h5ad', forbid_full_read)
    points, genes = result_metadata(path)
    assert points.cell_id.tolist() == data.obs_names.tolist()
    np.testing.assert_array_equal(points[['UMAP1', 'UMAP2']].to_numpy(), data.obsm['X_umap'])
    vector = gene_expression(path, 'id1')
    np.testing.assert_array_equal(vector, values[:, 1])
    assert find_genes(genes, 'genea').index.tolist() == ['id0', 'id1']
    assert find_genes(genes, 'id2').index.tolist() == ['id2']
    assert find_genes(genes, '').empty
    assert find_genes(genes, 'absent').empty
    summary = expression_summary(points, vector).set_index('cluster')
    assert summary.detected_percent.to_dict() == {'0': 50., '1': 50.}
    assert summary.mean_log1p_expression.to_dict() == {'0': 1., '1': 2.5}
    with pytest.raises(InputError, match='GENE_MISSING'):
        gene_expression(path, 'missing')
