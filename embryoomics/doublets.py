"""Optional raw-count Scrublet scoring, separated from human filtering decisions."""

import hashlib
import json
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .config import InputError


def score_doublets(data, settings, seed):
    # A fresh object prevents input X/raw or existing analysis metadata influencing scoring.
    if data.n_obs < 10 or data.n_vars < 3:
        raise InputError('SCRUBLET_SIZE: at least 10 QC-retained cells and 3 genes are required; disable Scrublet for tiny data')
    if 'batch' in data.obs and data.obs['batch'].nunique() > 1:
        raise InputError('SCRUBLET_BATCH: score one capture batch per run; split batches explicitly')
    work = ad.AnnData(data.layers['counts'].copy(), obs=data.obs[['sample_id']].copy())
    if (np.asarray((work.X > 0).sum(axis=1)).ravel() < 3).any():
        raise InputError('SCRUBLET_GENES: some cells have fewer than 3 detected genes; revise QC or disable Scrublet')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            sc.pp.scrublet(work, expected_doublet_rate=settings['expected_doublet_rate'],
                          threshold=settings['threshold'], n_prin_comps=settings['n_prin_comps'],
                          random_state=seed, use_approx_neighbors=False, mean_center=False,
                          sim_doublet_ratio=2.0, verbose=False)
        except (ValueError, RuntimeError, ImportError, ZeroDivisionError) as exc:
            raise InputError(f'SCRUBLET_FAILED: {exc}; check QC and reduce n_prin_comps if needed, or disable Scrublet') from exc
    values = work.obs['doublet_score'].to_numpy(dtype=float)
    simulated = np.asarray(work.uns['scrublet']['doublet_scores_sim'], dtype=float)
    if not np.isfinite(values).all() or not np.isfinite(simulated).all():
        raise InputError('SCRUBLET_UNSCORED: internal gene/cell filtering left unscored cells; revise QC or disable Scrublet')
    scores = pd.DataFrame({'doublet_score': values, 'predicted_doublet': values > settings['threshold']}, index=data.obs_names)
    scores.index.name = 'cell_id'
    summary = {'backend': 'scanpy.pp.scrublet', 'input_layer': 'counts (raw integer counts after ordinary QC)',
               'settings': settings.copy(), 'seed': seed, 'mean_center': False, 'normalize_variance': True,
               'log_transform': False, 'stdev_doublet_rate': 0.02, 'synthetic_doublet_umi_subsampling': 1.0,
               'knn_dist_metric': 'euclidean', 'n_neighbors': 'Scanpy default: round(0.5 * sqrt(n_obs))',
               'use_approx_neighbors': False, 'sim_doublet_ratio': 2.0,
               'scored_cells': len(scores), 'candidate_doublets': int(scores.predicted_doublet.sum()),
               'removed_cells': int(scores.predicted_doublet.sum()) if settings['mode'] == 'filter' else 0,
               'warnings': list(dict.fromkeys(str(w.message) for w in caught))}
    digest = hashlib.sha256(scores.to_csv(float_format='%.17g').encode())
    digest.update(json.dumps({k: v for k, v in summary.items() if k != 'warnings'}, sort_keys=True).encode())
    summary['score_digest'] = digest.hexdigest()
    return scores, simulated, summary


def apply_doublet_decision(metrics, qc_mask, scores, settings):
    """Keep all input cells in the audit; QC-failed cells have no doublet score."""
    audit = pd.DataFrame(index=metrics.index)
    audit.index.name = 'cell_id'
    audit['ordinary_qc_pass'] = qc_mask
    audit = audit.join(scores)
    audit['removed_as_doublet'] = False
    if settings['mode'] == 'filter':
        audit.loc[scores.index, 'removed_as_doublet'] = scores.predicted_doublet.to_numpy(dtype=bool)
    final = np.asarray(qc_mask) & ~audit.removed_as_doublet.to_numpy()
    audit['retained'] = final
    return final, audit
