"""Read result metadata and one expression column without loading the full matrix."""

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from .config import InputError


_SPARSE_CHUNK_SIZE = 1_000_000


def result_metadata(path):
    with h5py.File(path, 'r') as handle:
        obs = read_elem(handle['obs'])
        var = read_elem(handle['var'])
        coordinates = read_elem(handle['obsm']['X_umap'])
    if coordinates.shape != (len(obs), 2) or 'cluster' not in obs:
        raise InputError('RESULT_VIEW: expected two UMAP coordinates and cluster labels')
    points = obs[['cluster', 'sample_id', 'stage']].copy()
    points['cluster'] = points.cluster.astype(str)
    points['cell_id'] = points.index.astype(str)
    points['UMAP1'], points['UMAP2'] = coordinates[:, 0], coordinates[:, 1]
    genes = pd.DataFrame({'gene_id': var.index.astype(str),
                          'gene_symbol': var.get('gene_symbol', pd.Series('', index=var.index)).astype(str).to_numpy()})
    return points, genes.set_index('gene_id')


def find_genes(genes, query, limit=50):
    query = query.strip().lower()
    if not query:
        return genes.iloc[:0]
    ids = genes.index.str.lower()
    symbols = genes.gene_symbol.str.lower()
    selected = genes.loc[ids.str.contains(query, regex=False) | symbols.str.contains(query, regex=False)]
    exact = (selected.index.str.lower() == query) | (selected.gene_symbol.str.lower() == query)
    return pd.concat([selected.loc[exact], selected.loc[~exact]]).head(limit)


def gene_expression(path, gene_id):
    with h5py.File(path, 'r') as handle:
        var = read_elem(handle['var'])
        if gene_id not in var.index:
            raise InputError(f'GENE_MISSING: {gene_id}')
        column = int(var.index.get_loc(gene_id))
        storage = handle['X']
        if isinstance(storage, h5py.Dataset):
            return np.asarray(storage[:, column], dtype=float)
        encoding = storage.attrs.get('encoding-type')
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        shape = storage.attrs['shape']
        values = np.zeros(int(shape[0]), dtype=float)
        if encoding == 'csc_matrix':
            start, end = storage['indptr'][column:column + 2]
            np.add.at(values, storage['indices'][start:end], storage['data'][start:end])
        elif encoding == 'csr_matrix':
            pointers = storage['indptr'][()]
            # Bound temporary memory even when the saved result has many nonzero entries.
            for start in range(0, len(storage['indices']), _SPARSE_CHUNK_SIZE):
                indices = storage['indices'][start:start + _SPARSE_CHUNK_SIZE]
                locations = np.flatnonzero(indices == column) + start
                if len(locations):
                    rows = np.searchsorted(pointers[1:], locations, side='right')
                    np.add.at(values, rows, storage['data'][locations])
        else:
            raise InputError('RESULT_VIEW: unsupported expression matrix encoding')
        return values


def expression_summary(points, values):
    if len(points) != len(values):
        raise InputError('RESULT_VIEW: expression and cell metadata dimensions differ')
    summary = pd.DataFrame({'cluster': points.cluster.to_numpy(), 'expression': values, 'detected': values > 0})
    return summary.groupby('cluster', observed=True).agg(
        cells=('expression', 'size'), mean_log1p_expression=('expression', 'mean'),
        detected_percent=('detected', lambda series: 100 * series.mean())).reset_index()
