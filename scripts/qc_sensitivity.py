"""Run two explicit QC sensitivity scenarios without changing the baseline."""

import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml
from sklearn.metrics import adjusted_rand_score

from embryoomics.pipeline import run_analysis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    baseline = Path(args.baseline)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    original = ad.read_h5ad(baseline / 'processed.h5ad', backed='r')
    try:
        labels = original.obs['cluster'].astype(str).copy()
    finally:
        original.file.close()
    base_config = yaml.safe_load((baseline / 'config.yaml').read_text(encoding='utf-8'))
    rows = []
    for name, min_genes, max_mt in [('min200_mt20', 200, 20), ('min500_mt10', 500, 10)]:
        config = dict(base_config)
        config['qc'] = {'min_genes': min_genes, 'max_genes': None, 'max_pct_mt': max_mt}
        config_path = out / f'{name}.yaml'
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding='utf-8')
        result = run_analysis(config_path, args.input, out / name)
        processed = ad.read_h5ad(out / name / 'processed.h5ad', backed='r')
        try:
            new = processed.obs['cluster'].astype(str).copy()
        finally:
            processed.file.close()
        shared = labels.loc[new.index]
        contingency = pd.crosstab(shared, new)
        contingency.to_csv(out / f'{name}_cluster_overlap.csv')
        retention = pd.DataFrame({'baseline_cells': labels.value_counts(), 'retained_cells': shared.value_counts()}).fillna(0).astype(int)
        retention['retained_percent'] = 100 * retention.retained_cells / retention.baseline_cells
        retention.to_csv(out / f'{name}_baseline_cluster_retention.csv', index_label='baseline_cluster')
        row = {**result, 'scenario': name, 'qc': config['qc'],
               'retained_percent': round(100 * len(new) / len(labels), 2),
               'ari_on_shared_cells': adjusted_rand_score(shared, new),
               'interpretation': 'Exploratory QC sensitivity; ARI compares shared cells only and is not biological accuracy.'}
        rows.append(row)
        print(json.dumps(row), flush=True)
    (out / 'comparison.json').write_text(json.dumps(rows, indent=2) + '\n', encoding='utf-8')
    markers = pd.read_csv(baseline / 'markers.csv', dtype={'cluster_id': str})
    markers.loc[markers.cluster_id.isin(['2', '3', '7'])].to_csv(out / 'baseline_flagged_cluster_markers.csv', index=False)


if __name__ == '__main__':
    main()
