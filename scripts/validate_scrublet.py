"""Software sensitivity experiment, not a recommendation for biological filtering.

Run from the repository root with the installed development package.
"""
import json
import argparse
from pathlib import Path
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
import yaml

from embryoomics.pipeline import inspect_doublets, run_analysis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='results/GSM8559287_scrublet_acceptance_20261113')
    root = Path(parser.parse_args().out)
    root.mkdir(parents=True, exist_ok=False)
    source = Path('examples/GSM8559287.yaml')
    config = yaml.safe_load(source.read_text(encoding='utf-8'))
    config['doublets'] = {'mode': 'score', 'expected_doublet_rate': 0.05, 'threshold': 0.15,
                         'n_prin_comps': 30, 'upstream_handling': 'unknown'}
    path = root / 'score.yaml'
    path.write_text(yaml.safe_dump(config), encoding='utf-8')
    preview = inspect_doublets('data/GSM8559287', path)
    preview['audit'].to_csv(root / 'preview_scores.csv')
    run_analysis(path, 'data/GSM8559287', root / 'score', expected_doublet_digest=preview['summary']['score_digest'])
    config['doublets']['mode'] = 'filter'
    path = root / 'filter.yaml'
    path.write_text(yaml.safe_dump(config), encoding='utf-8')
    filtered_preview = inspect_doublets('data/GSM8559287', path)
    run_analysis(path, 'data/GSM8559287', root / 'filter', confirm_doublets=True,
                 expected_doublet_digest=filtered_preview['summary']['score_digest'])
    score = ad.read_h5ad(root / 'score/processed.h5ad')
    filtered = ad.read_h5ad(root / 'filter/processed.h5ad')
    baseline = ad.read_h5ad('results/GSM8559287_scanpy_plots_20261113/processed.h5ad')
    assert score.obs_names.equals(baseline.obs_names)
    assert (score.layers['counts'] != baseline.layers['counts']).nnz == 0
    assert (filtered.layers['counts'] != baseline[filtered.obs_names].layers['counts']).nnz == 0
    score_audit = pd.read_csv(root / 'score/doublet_scores.csv', index_col='cell_id')
    filter_audit = pd.read_csv(root / 'filter/doublet_scores.csv', index_col='cell_id')
    np.testing.assert_allclose(score_audit.doublet_score, filter_audit.doublet_score)
    assert score.n_obs == len(score_audit) and not score_audit.removed_as_doublet.any()
    assert filtered.n_obs == int(filter_audit.retained.sum()) == filtered_preview['retained']
    for mode in ('score', 'filter'):
        manifest = json.loads((root / mode / 'run_manifest.json').read_text(encoding='utf-8'))
        with zipfile.ZipFile(root / mode / 'report_bundle.zip') as bundle:
            assert set(bundle.namelist()) == set(manifest['files']) and bundle.testzip() is None
    markers_before = pd.read_csv(root / 'score/markers.csv').groupby('cluster_id').head(10)
    markers_after = pd.read_csv(root / 'filter/markers.csv').groupby('cluster_id').head(10)
    # Cluster IDs can change; ARI compares partitions and marker overlap is only descriptive.
    summary = {'purpose': 'software sensitivity only; rate 0.05 and threshold 0.15 are not scientific recommendations',
               'score_cells': score.n_obs, 'filtered_cells': filtered.n_obs,
               'candidates': int(filter_audit.predicted_doublet.sum()),
               'score_vs_baseline_ari': adjusted_rand_score(baseline.obs.cluster, score.obs.cluster),
               'filter_vs_score_common_cells_ari': adjusted_rand_score(score[filtered.obs_names].obs.cluster, filtered.obs.cluster),
               'score_clusters': score.obs.cluster.nunique(), 'filtered_clusters': filtered.obs.cluster.nunique(),
               'top10_marker_unique_gene_overlap': len(set(markers_before.gene_id) & set(markers_after.gene_id)),
               'score_top10_unique_genes': markers_before.gene_id.nunique(),
               'filter_top10_unique_genes': markers_after.gene_id.nunique(),
               'preview_matches_run': True, 'counts_preserved': True, 'zip_integrity': True}
    (root / 'acceptance.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
