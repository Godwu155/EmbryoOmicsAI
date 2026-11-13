"""Summarize observed marker expression across QC scenarios, without annotation."""

import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


PANEL = ['Ttr', 'Mixl1', 'Sox17', 'Foxa2', 'Lhx1', 'Sp5', 'T', 'Pou5f1', 'Dnmt3a', 'Dnmt3b', 'Spink1', 'Apoe']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', nargs='+', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    missing = {}
    fig, axes = plt.subplots(len(args.runs), 1, figsize=(13, 4 * len(args.runs)), squeeze=False)
    for run_path, axis in zip(args.runs, axes[:, 0]):
        run = Path(run_path)
        data = ad.read_h5ad(run / 'processed.h5ad')
        symbols = data.var.gene_symbol.astype(str).to_numpy()
        labels = data.obs.cluster.astype(str).to_numpy()
        clusters = sorted(set(labels), key=lambda x: int(x) if x.isdigit() else x)
        missing[run.name] = [gene for gene in PANEL if gene not in symbols]
        for gene_idx, gene in enumerate(PANEL):
            columns = np.flatnonzero(symbols == gene)
            if not len(columns):
                continue
            # Preserve ambiguity: aggregate counts for duplicate symbols and record ID count.
            totals = np.asarray(data.layers['counts'][:, columns].sum(axis=1)).ravel()
            for cluster_idx, cluster in enumerate(clusters):
                selected = labels == cluster
                fraction = float((totals[selected] > 0).mean())
                rows.append({'run': run.name, 'cluster': cluster, 'gene_symbol': gene,
                             'matching_gene_ids': len(columns), 'cells': int(selected.sum()),
                             'fraction_detected': fraction, 'mean_raw_counts': float(totals[selected].mean())})
                axis.scatter(gene_idx, cluster_idx, s=200 * fraction, c=[fraction], cmap='viridis', vmin=0, vmax=1)
        axis.set(xticks=range(len(PANEL)), xticklabels=PANEL, yticks=range(len(clusters)),
                 yticklabels=clusters, ylabel='Run-local cluster', title=run.name)
        axis.tick_params(axis='x', rotation=45)
    fig.suptitle('Observed marker detection — dot area and color represent fraction detected; no cell-type assignments')
    fig.tight_layout(rect=(0, 0, 1, .96))
    fig.savefig(out / 'marker_detection.png', dpi=150)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out / 'marker_detection.csv', index=False)
    pd.Series({key: ','.join(value) for key, value in missing.items()}, name='missing_symbols').to_csv(out / 'missing_symbols.csv')


if __name__ == '__main__':
    main()
