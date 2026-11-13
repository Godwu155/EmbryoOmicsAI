"""Read-only live API acceptance; run with the isolated CELLxGENE interpreter."""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import requests
from server.common.fbs.matrix import decode_matrix_fbs


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exchange', required=True)
    parser.add_argument('--url', default='http://127.0.0.1:5057')
    args = parser.parse_args()
    root = Path(args.exchange)
    provenance = json.loads((root / 'exchange.json').read_text(encoding='utf-8'))
    api = args.url.rstrip('/') + '/api/v0.2'
    def get(path, params=None, binary=False):
        response = requests.get(api + path, params=params,
            headers={'Accept': 'application/octet-stream'} if binary else {}, timeout=30)
        response.raise_for_status()
        return decode_matrix_fbs(response.content) if binary else response.json()
    schema = get('/schema')['schema']
    obs = get('/annotations/obs', [('annotation-name', 'cluster'), ('annotation-name', 'candidate_cell_type')], True)
    var = get('/annotations/var', [('annotation-name', 'viewer_gene'), ('annotation-name', 'gene_id')], True)
    layout = get('/layout/obs', {'layout-name': 'umap'}, True)
    view = ad.read_h5ad(root / 'viewer.h5ad')
    source = ad.read_h5ad(Path(provenance['source_run']) / 'processed.h5ad', backed='r')
    try:
        assert len(obs) == view.n_obs and len(var) == view.n_vars
        assert obs.cluster.astype(str).tolist() == view.obs.cluster.astype(str).tolist()
        assert var.viewer_gene.astype(str).tolist() == view.var.viewer_gene.astype(str).tolist()
        assert var.viewer_gene.is_unique
        assert layout.shape == (view.n_obs, 2) and np.isfinite(layout.to_numpy()).all()
        np.testing.assert_array_equal(view.obsm['X_umap'], source.obsm['X_umap'])
        gene = 'ENSMUSG00000061808'
        label = view.var.loc[gene, 'viewer_gene']
        expression = get('/data/var', {'var:viewer_gene': label}, True).to_numpy().ravel()
        original = np.asarray(source.X[:, source.var_names.get_loc(gene)].toarray()).ravel()
        np.testing.assert_allclose(expression, original, rtol=1e-6, atol=1e-6)
        assert digest(Path(provenance['source_run']) / 'processed.h5ad') == provenance['processed_sha256']
        assert digest(Path(provenance['source_run']) / 'annotations.csv') == provenance['annotations_sha256']
        summary = {'cells': len(obs), 'genes': len(var), 'clusters': obs.cluster.nunique(),
            'umap_dimensions': list(layout.shape), 'embedding_preserved': True,
            'gene_tested': label, 'expression_max_abs_error': float(np.max(np.abs(expression - original))),
            'unique_display_labels': True, 'duplicate_symbol_genes': int(view.var.gene_symbol.duplicated(keep=False).sum()),
            'source_unchanged': True, 'viewer_bytes': (root / 'viewer.h5ad').stat().st_size,
            'source_bytes': (Path(provenance['source_run']) / 'processed.h5ad').stat().st_size,
            'editable_columns': [column['name'] for column in schema['annotations']['obs']['columns'] if column.get('writable')]}
        (root / 'api_acceptance.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
        print(json.dumps(summary, indent=2))
    finally:
        source.file.close()


if __name__ == '__main__':
    main()
