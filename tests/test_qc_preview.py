import numpy as np
import pandas as pd
import pytest

from embryoomics.config import InputError
from embryoomics.qc import preview_qc_metrics


def test_qc_rule_overlap_and_sample_totals():
    metrics = pd.DataFrame({'sample_id': ['a'] * 3 + ['b'] * 3,
                            'n_genes_by_counts': [0, 2, 4, 6, 8, 10],
                            'pct_counts_mt': [0, 50, 30, 10, np.nan, 5]})
    mask, summary, reasons = preview_qc_metrics(metrics, {'min_genes': 3, 'max_genes': 9, 'max_pct_mt': 20})
    assert mask.tolist() == [False, False, False, True, False, False]
    assert summary.set_index('sample_id').retained_cells.to_dict() == {'a': 0, 'b': 1}
    assert reasons.failed_cells.tolist() == [2, 1, 3]
    assert reasons.newly_removed_cells.tolist() == [2, 1, 2]
    assert reasons.retained_after_rule.tolist() == [4, 3, 1]
    assert reasons.newly_removed_cells.sum() == summary.removed_cells.sum() == 5


def test_disabled_qc_rules_and_missing_mitochondrial_values():
    metrics = pd.DataFrame({'sample_id': ['a', 'a'], 'n_genes_by_counts': [2, 4], 'pct_counts_mt': [np.nan, np.nan]})
    mask, _, reasons = preview_qc_metrics(metrics, {'min_genes': None, 'max_genes': None, 'max_pct_mt': None})
    assert mask.all() and not reasons.enabled.any()
    with pytest.raises(InputError, match='MT_UNAVAILABLE'):
        preview_qc_metrics(metrics, {'min_genes': None, 'max_genes': None, 'max_pct_mt': 20})
