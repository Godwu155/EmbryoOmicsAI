"""Cluster marker evidence and human review records."""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from .config import InputError


ANNOTATION_COLUMNS = ["cluster_id", "candidate_cell_type", "support_genes", "conflict_genes", "reference_url", "review_status", "reviewer_note"]


def marker_table(data, top_n=30):
    clusters = sorted(data.obs["cluster"].astype(str).unique())
    if len(clusters) < 2:
        return pd.DataFrame(columns=["cluster_id", "gene_id", "gene_symbol", "pct_in", "pct_reference", "logfoldchange", "pval_adj", "test_method"])
    sc.tl.rank_genes_groups(data, groupby="cluster", method="wilcoxon", use_raw=False, n_genes=min(top_n, data.n_vars), pts=True)
    ranking = sc.get.rank_genes_groups_df(data, group=None)
    rows = []
    for _, row in ranking.iterrows():
        cluster = str(row["group"])
        symbol = str(row["names"])
        gene_idx = data.var_names.get_loc(symbol)
        inside = data.obs["cluster"].astype(str).to_numpy() == cluster
        raw = data.layers["counts"][:, gene_idx]
        present = np.asarray((raw > 0).sum(axis=1)).ravel() if sparse.issparse(raw) else np.asarray(raw > 0).ravel()
        rows.append({"cluster_id": cluster, "gene_id": symbol, "gene_symbol": str(data.var.iloc[gene_idx]["gene_symbol"]),
                     "pct_in": float(present[inside].mean()), "pct_reference": float(present[~inside].mean()),
                     "logfoldchange": float(row["logfoldchanges"]), "pval_adj": float(row["pvals_adj"]), "test_method": "Scanpy Wilcoxon on log1p normalized expression; cluster vs rest"})
    return pd.DataFrame(rows)


def load_annotations(path, clusters):
    if path is None:
        return pd.DataFrame([{**{column: "" for column in ANNOTATION_COLUMNS}, "cluster_id": cluster, "candidate_cell_type": "unknown", "review_status": "unreviewed"} for cluster in clusters])
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = set(ANNOTATION_COLUMNS) - set(table)
    if missing:
        raise InputError(f"ANNOTATIONS: missing columns {sorted(missing)}")
    if table.cluster_id.duplicated().any() or set(table.cluster_id) != set(clusters):
        raise InputError("ANNOTATIONS: exactly one row per current cluster is required")
    if (table.candidate_cell_type == "").any():
        raise InputError("ANNOTATIONS: use 'unknown' for unresolved cell types")
    return table[ANNOTATION_COLUMNS].sort_values("cluster_id")
