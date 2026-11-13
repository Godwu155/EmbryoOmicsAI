import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from embryoomics.reporting import save_marker_figures


def test_duplicate_symbols_and_raw_layer_do_not_change_plot_statistics(tmp_path, monkeypatch):
    data = ad.AnnData(sparse.csr_matrix([[0., 2.], [0., 4.], [1., 0.], [3., 0.]]))
    data.var_names = ['id1', 'id2']
    data.obs['cluster'] = ['0', '0', '1', '1']
    data.layers['counts'] = sparse.csr_matrix(np.ones((4, 2)) * 100)
    data.raw = data.copy()
    markers = pd.DataFrame({'cluster_id': ['0', '1'], 'gene_id': ['id2', 'id1'],
                            'gene_symbol': ['Duplicate', 'Duplicate'], 'logfoldchange': [2., 3.]})
    original = data.X.copy()
    dotplot = sc.pl.dotplot
    def checked_dotplot(*args, **kwargs):
        plot = dotplot(*args, **kwargs)
        np.testing.assert_allclose(plot.dot_color_df.to_numpy(), [[3., 0.], [0., 2.]])
        np.testing.assert_allclose(plot.dot_size_df.to_numpy(), [[1., 0.], [0., 1.]])
        assert plot.dot_color_df.columns.tolist() == ['Duplicate (id2)', 'Duplicate (id1)']
        return plot
    monkeypatch.setattr(sc.pl, 'dotplot', checked_dotplot)
    details = save_marker_figures(data, tmp_path, markers)
    assert details['gene_ids'] == ['id2', 'id1']
    assert details['status'] == 'generated'
    assert all((tmp_path / name).is_file() for name in details['files'])
    np.testing.assert_array_equal(data.X.toarray(), original.toarray())
    assert data.var_names.tolist() == ['id1', 'id2']


def test_no_positive_markers_skips_figures(tmp_path):
    data = ad.AnnData(sparse.csr_matrix([[1.], [2.]]))
    markers = pd.DataFrame(columns=['cluster_id', 'gene_id', 'gene_symbol', 'logfoldchange'])
    details = save_marker_figures(data, tmp_path, markers)
    assert details['status'] == 'skipped_no_positive_markers'
    assert details['files'] == ['marker_plot_genes.csv']
    assert pd.read_csv(tmp_path / 'marker_plot_genes.csv').empty
