import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from embryoomics.qc import calculate_qc
from embryoomics.config import InputError


def test_explicit_gene_ids_override_symbol_guess():
    data = ad.AnnData(sparse.csr_matrix([[1, 3], [4, 0]]))
    data.var_names = ['a', 'b']
    data.var['gene_symbol'] = pd.Categorical(['mt-Fallback', 'Other'])
    data.layers['counts'] = data.X.copy()
    method, _ = calculate_qc(data, ['b'])
    assert method['method'] == 'configured mitochondrial gene IDs'
    np.testing.assert_array_equal(data.obs.pct_counts_mt, [75., 0.])
    calculate_qc(data)
    np.testing.assert_array_equal(data.obs.pct_counts_mt, [25., 100.])
    with pytest.raises(InputError, match='MT_GENE_IDS'):
        calculate_qc(data, ['missing'])
