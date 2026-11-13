import json
import gzip
import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import hypergeom, false_discovery_control
import yaml

from embryoomics.enrichment import enrich_markers
from embryoomics.config import InputError


@pytest.fixture
def enrichment_input(tmp_path):
    run = tmp_path / 'run'
    run.mkdir()
    data = ad.AnnData(sparse.csr_matrix(np.ones((4, 11))))
    data.var_names = [f'g{i}' for i in range(11)]
    data.var['gene_symbol'] = [f'S{i}' for i in range(8)] + ['Same', 'Same', 'Absent']
    data.obs['cluster'] = ['0', '0', '1', '1']
    counts = np.ones((4, 11), dtype=int)
    counts[:, 10] = 0
    data.layers['counts'] = sparse.csr_matrix(counts)
    data.write_h5ad(run / 'processed.h5ad')
    (run / 'run_manifest.json').write_text(json.dumps({'species': 'mus_musculus'}), encoding='utf-8')
    pd.DataFrame({'cluster_id': ['0', '0', '0', '1'], 'gene_id': ['g0', 'g1', 'g8', 'g4'],
        'logfoldchange': [2., 3., 2., 4.], 'pval_adj': [.001] * 4}).to_csv(run / 'markers.csv', index=False)
    gmt = tmp_path / 'sets.gmt'
    gmt.write_text('A\tdesc\tS0\tS1\nB\tdesc\tS6\tS7\nC\tdesc\tS3\tS4\n', encoding='utf-8')
    meta = tmp_path / 'metadata.yaml'
    meta.write_text(yaml.safe_dump({'species': 'mus_musculus', 'id_type': 'gene_symbol',
        'source_url': 'https://example.invalid/synthetic', 'version': 'synthetic', 'license': 'synthetic fixture'}), encoding='utf-8')
    return run, gmt, meta


def test_offline_background_mapping_and_full_fdr_family(enrichment_input, tmp_path, monkeypatch):
    import requests
    def forbid_network(*args, **kwargs):
        pytest.fail('Enrichment must remain offline')
    monkeypatch.setattr(requests.sessions.Session, 'request', forbid_network)
    run, gmt, meta = enrichment_input
    before = (run / 'processed.h5ad').read_bytes()
    out = tmp_path / 'enrichment'
    summary = enrich_markers(run, gmt, meta, out, min_size=2)
    assert summary['background_size'] == 8 and summary['cluster_term_tests'] == 6
    results = pd.read_csv(out / 'enrichment.csv', dtype={'cluster_id': str}).set_index(['cluster_id', 'term'])
    assert results.loc[('0', 'A'), 'pvalue'] == pytest.approx(hypergeom.sf(1, 8, 2, 2))
    assert results.loc[('0', 'B'), 'pvalue'] == 1
    np.testing.assert_allclose(results.pval_adj, false_discovery_control(results.pvalue.to_numpy()))
    mapping = pd.read_csv(out / 'gene_mapping.csv')
    assert not mapping.loc[mapping.gene_id.isin(['g8', 'g9', 'g10']), 'included'].any()
    assert (run / 'processed.h5ad').read_bytes() == before


def test_wrong_species_and_changed_snapshot_rejected(enrichment_input, tmp_path):
    run, gmt, meta = enrichment_input
    config = yaml.safe_load(meta.read_text())
    config['species'] = 'homo_sapiens'
    meta.write_text(yaml.safe_dump(config))
    with pytest.raises(InputError, match='SPECIES'):
        enrich_markers(run, gmt, meta, tmp_path / 'wrong')
    config.update(species='mus_musculus', gmt_sha256='changed')
    meta.write_text(yaml.safe_dump(config))
    with pytest.raises(InputError, match='GENE_SET_HASH'):
        enrich_markers(run, gmt, meta, tmp_path / 'changed')


def test_gaf_preparation_excludes_negative_other_taxa_and_ambiguous_symbols(tmp_path):
    spec = importlib.util.spec_from_file_location('prepare_go', Path(__file__).resolve().parents[1] / 'scripts/prepare_mouse_go.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    path = tmp_path / 'test.gaf.gz'
    def row(symbol, gene, qualifier='', taxon='NCBITaxon:10090'):
        fields = ['MGI', gene, symbol, qualifier, 'GO:1', 'PMID:1', 'IMP', '', 'P', '', '', 'protein', taxon, '', 'MGI', '', '']
        return '\t'.join(fields) + '\n'
    with gzip.open(path, 'wt') as handle:
        handle.write('!date-generated: test\n' + row('Good', '1') + row('Negative', '2', 'NOT') +
            row('Other', '3', taxon='taxon:9606') + row('Ambiguous', '4') + row('Ambiguous', '5'))
    module.prepare(path, tmp_path / 'sets')
    assert (tmp_path / 'sets/mouse_GO_BP_direct.gmt').read_text().strip().endswith('\tGood')
