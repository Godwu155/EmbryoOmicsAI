"""QC metrics and explicit filtering preview."""

import numpy as np
import pandas as pd

from .config import InputError


def calculate_qc(data, mitochondrial_gene_ids=None):
    counts = data.layers["counts"].tocsr()
    data.obs["total_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    data.obs["n_genes_by_counts"] = np.diff(counts.indptr)
    symbols = data.var["gene_symbol"].astype('string').fillna("").astype(str)
    is_mt = data.var_names.isin(mitochondrial_gene_ids) if mitochondrial_gene_ids is not None else symbols.str.match(r"^mt-", case=False).to_numpy()
    if mitochondrial_gene_ids is not None and not set(mitochondrial_gene_ids).issubset(set(data.var_names)):
        raise InputError('MT_GENE_IDS: configured mitochondrial gene IDs are missing from the input')
    warnings = []
    if is_mt.any():
        mt_counts = np.asarray(counts[:, is_mt].sum(axis=1)).ravel()
        totals = data.obs["total_counts"].to_numpy()
        data.obs["pct_counts_mt"] = np.divide(mt_counts * 100, totals, out=np.full(data.n_obs, np.nan), where=totals > 0)
        method = 'configured mitochondrial gene IDs' if mitochondrial_gene_ids is not None else "provided gene symbols, case-insensitive ^mt-"
    else:
        data.obs["pct_counts_mt"] = np.nan
        method = "no matching gene symbols"
        warnings.append("No mitochondrial gene symbols matched; pct_counts_mt unavailable and mitochondrial filtering disabled.")
    return {"method": method, "matched_genes": int(is_mt.sum())}, warnings


def preview_qc(data, thresholds):
    mask, summary, _ = preview_qc_metrics(data.obs, thresholds)
    return mask, summary


def preview_qc_metrics(metrics, thresholds):
    """Preview the exact filtering rules from cached metrics, without loading counts."""
    mask = np.ones(len(metrics), dtype=bool)
    genes = metrics['n_genes_by_counts'].to_numpy()
    mt = metrics['pct_counts_mt'].to_numpy()
    reasons = []
    for key in ('min_genes', 'max_genes', 'max_pct_mt'):
        threshold = thresholds[key]
        passed = np.ones(len(metrics), dtype=bool)
        if threshold is not None:
            if key == 'min_genes':
                passed = genes >= threshold
            elif key == 'max_genes':
                passed = genes <= threshold
            else:
                if not np.isfinite(mt).any():
                    raise InputError('MT_UNAVAILABLE: no mitochondrial genes matched; set qc.max_pct_mt: null')
                passed = np.isfinite(mt) & (mt <= threshold)
        newly_removed = int((mask & ~passed).sum())
        mask &= passed
        reasons.append({'filter_rule': key, 'threshold': threshold, 'enabled': threshold is not None,
                        'failed_cells': int((~passed).sum()), 'newly_removed_cells': newly_removed,
                        'retained_after_rule': int(mask.sum())})
    preview = pd.DataFrame({"sample_id": metrics["sample_id"].astype(str), "retained": mask})
    summary = preview.groupby("sample_id", sort=True).agg(input_cells=("retained", "size"), retained_cells=("retained", "sum")).reset_index()
    summary["removed_cells"] = summary["input_cells"] - summary["retained_cells"]
    return mask, summary, pd.DataFrame(reasons)
