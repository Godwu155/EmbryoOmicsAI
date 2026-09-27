"""QC metrics and explicit filtering preview."""

import numpy as np
import pandas as pd

from .config import InputError


def calculate_qc(data):
    counts = data.layers["counts"].tocsr()
    data.obs["total_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    data.obs["n_genes_by_counts"] = np.diff(counts.indptr)
    symbols = data.var["gene_symbol"].fillna("").astype(str)
    is_mt = symbols.str.match(r"^mt-", case=False).to_numpy()
    warnings = []
    if is_mt.any():
        mt_counts = np.asarray(counts[:, is_mt].sum(axis=1)).ravel()
        totals = data.obs["total_counts"].to_numpy()
        data.obs["pct_counts_mt"] = np.divide(mt_counts * 100, totals, out=np.full(data.n_obs, np.nan), where=totals > 0)
        method = "mouse gene symbols, case-insensitive ^mt-"
    else:
        data.obs["pct_counts_mt"] = np.nan
        method = "no matching mouse gene symbols"
        warnings.append("No mitochondrial gene symbols matched; pct_counts_mt unavailable and mitochondrial filtering disabled.")
    return {"method": method, "matched_genes": int(is_mt.sum())}, warnings


def preview_qc(data, thresholds):
    mask = np.ones(data.n_obs, dtype=bool)
    genes = data.obs["n_genes_by_counts"].to_numpy()
    if thresholds["min_genes"] is not None:
        mask &= genes >= thresholds["min_genes"]
    if thresholds["max_genes"] is not None:
        mask &= genes <= thresholds["max_genes"]
    if thresholds["max_pct_mt"] is not None:
        values = data.obs["pct_counts_mt"].to_numpy()
        if not np.isfinite(values).any():
            raise InputError("MT_UNAVAILABLE: no mitochondrial genes matched; set qc.max_pct_mt: null")
        mask &= np.isfinite(values) & (values <= thresholds["max_pct_mt"])
    preview = pd.DataFrame({"sample_id": data.obs["sample_id"].astype(str), "retained": mask})
    summary = preview.groupby("sample_id", sort=True).agg(input_cells=("retained", "size"), retained_cells=("retained", "sum")).reset_index()
    summary["removed_cells"] = summary["input_cells"] - summary["retained_cells"]
    return mask, summary
