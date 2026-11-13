"""Recompute input QC and compare it with a saved run; never change thresholds."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad

from embryoomics.io import checksums
from embryoomics.pipeline import inspect_input


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    out = Path(args.out)
    if out.exists():
        raise ValueError('Choose a new audit output directory')
    baseline = Path(args.baseline)
    manifest = json.loads((baseline / 'run_manifest.json').read_text(encoding='utf-8'))
    info = inspect_input(args.input, args.config)
    metrics = info['qc_metrics']
    fresh_hashes = {Path(item['file']).name: item['sha256'] for item in checksums(info['files'])}
    saved_hashes = {Path(item['file']).name: item['sha256'] for item in manifest['input_files']}
    assert fresh_hashes == saved_hashes, 'Input checksums differ from baseline'
    saved = pd.read_csv(baseline / 'qc_metrics.csv', index_col='cell_id')
    pd.testing.assert_frame_equal(metrics, saved, check_dtype=False, check_names=False, rtol=1e-10, atol=1e-10)
    missing = [name for name in manifest['files'] if not (baseline / name).is_file()]
    assert not missing, f'Missing baseline artifacts: {missing}'
    quantiles = metrics[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].quantile([0, .01, .05, .25, .5, .75, .95, .99, 1])
    rows = []
    # These values explore sensitivity, not biologically approved thresholds.
    for min_genes in (None, 200, 500, 1000):
        for max_mt in (None, 5, 10, 20):
            mask = pd.Series(True, index=metrics.index)
            if min_genes is not None:
                mask &= metrics.n_genes_by_counts >= min_genes
            if max_mt is not None:
                mask &= metrics.pct_counts_mt.notna() & (metrics.pct_counts_mt <= max_mt)
            rows.append({'min_genes': min_genes, 'max_pct_mt': max_mt, 'input_cells': len(metrics),
                         'retained_cells': int(mask.sum()), 'removed_cells': int((~mask).sum()),
                         'retained_percent': round(100 * mask.mean(), 2)})
    out.mkdir(parents=True)
    quantiles.to_csv(out / 'qc_quantiles.csv', index_label='quantile')
    pd.DataFrame(rows).to_csv(out / 'threshold_sensitivity.csv', index=False)
    processed = ad.read_h5ad(baseline / 'processed.h5ad', backed='r')
    try:
        cluster_qc = processed.obs.groupby('cluster', observed=True).agg(
            cells=('total_counts', 'size'), median_genes=('n_genes_by_counts', 'median'),
            median_mt=('pct_counts_mt', 'median'),
            high_mt_percent=('pct_counts_mt', lambda x: 100 * (x > 10).mean()),
            low_genes_percent=('n_genes_by_counts', lambda x: 100 * (x < 500).mean()))
        cluster_qc.to_csv(out / 'cluster_qc.csv')
    finally:
        processed.file.close()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].hist(metrics.n_genes_by_counts, bins=60)
    axes[0].set(xlabel='Detected genes per cell', ylabel='Cells')
    axes[1].hist(metrics.pct_counts_mt, bins=60)
    axes[1].set(xlabel='Mitochondrial counts (%)', ylabel='Cells')
    points = axes[2].scatter(metrics.total_counts, metrics.n_genes_by_counts, c=metrics.pct_counts_mt, s=3, alpha=.5)
    fig.colorbar(points, ax=axes[2], label='Mitochondrial counts (%)')
    axes[2].set(xscale='log', xlabel='Total counts (log scale)', ylabel='Detected genes')
    fig.suptitle('GSM8559287 input QC — no additional filtering applied')
    fig.tight_layout()
    fig.savefig(out / 'qc_review.png', dpi=150)
    plt.close(fig)
    result = {'sample_id': manifest['sample_id'], 'cells': info['cells'], 'genes': info['genes'],
              'input_hashes_match': True, 'recomputed_qc_matches': True, 'all_baseline_artifacts_present': True,
              'baseline_thresholds': manifest['parameters']['qc'], 'mt_detection': info['mt'],
              'review_required': 'Sensitivity grid is illustrative; thresholds and cell annotations require scientific review.',
              'input_sha256': fresh_hashes}
    (out / 'audit.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
