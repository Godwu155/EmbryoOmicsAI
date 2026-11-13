"""Offline marker over-representation using an explicit local gene-set snapshot."""
import argparse
import json
from pathlib import Path
from html import escape
import importlib.metadata
import zipfile

import anndata as ad
import gseapy as gp
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import yaml

from .cellxgene import sha256
from .config import InputError


def read_gmt(path):
    sets = {}
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 3 or not fields[0] or fields[0] in sets:
            raise InputError('GMT: unique named terms, a description and at least one gene are required')
        sets[fields[0]] = set(filter(None, fields[2:]))
    if not sets:
        raise InputError('GMT: no gene sets')
    return sets


def enrich_markers(run, gmt, metadata, out, *, marker_fdr=0.05, min_size=5, max_size=500):
    run, out = Path(run).resolve(), Path(out).resolve()
    if out.exists():
        raise InputError('OUTPUT_EXISTS: choose a new enrichment directory')
    if not 0 < marker_fdr <= 1 or not 1 <= min_size <= max_size:
        raise InputError('ENRICHMENT_CONFIG: invalid FDR or gene-set size range')
    meta = yaml.safe_load(Path(metadata).read_text(encoding='utf-8'))
    for key in ('species', 'id_type', 'source_url', 'version', 'license'):
        if not isinstance(meta, dict) or not isinstance(meta.get(key), str) or not meta[key].strip():
            raise InputError(f'GENE_SET_METADATA: explicit {key} is required')
    if meta['id_type'] not in ('gene_id', 'gene_symbol'):
        raise InputError('GENE_SET_METADATA: id_type must be gene_id or gene_symbol')
    manifest = json.loads((run / 'run_manifest.json').read_text(encoding='utf-8'))
    if meta['species'] != manifest['species']:
        raise InputError('GENE_SET_SPECIES: gene-set species must match the run')
    inputs = {name: sha256(path) for name, path in {'processed': run / 'processed.h5ad', 'markers': run / 'markers.csv',
              'manifest': run / 'run_manifest.json', 'gmt': gmt, 'metadata': metadata}.items()}
    if meta.get('gmt_sha256') and meta['gmt_sha256'] != inputs['gmt']:
        raise InputError('GENE_SET_HASH: GMT differs from its declared snapshot')
    sets = read_gmt(gmt)
    data = ad.read_h5ad(run / 'processed.h5ad', backed='r')
    try:
        if 'counts' not in data.layers:
            raise InputError('ENRICHMENT_COUNTS: processed data needs raw counts')
        detected = np.zeros(data.n_vars, dtype=bool)
        for start in range(0, data.n_obs, 256):
            detected |= np.asarray((data.layers['counts'][start:start + 256] > 0).sum(axis=0)).ravel() > 0
        mapping = pd.DataFrame({'gene_id': data.var_names.astype(str),
            'gene_symbol': data.var.gene_symbol.astype('string').fillna('').astype(str).to_numpy(), 'detected': detected})
        mapping['ambiguous_symbol'] = mapping.gene_symbol.duplicated(keep=False) | mapping.gene_symbol.eq('')
        mapping['included'] = detected & (~mapping.ambiguous_symbol if meta['id_type'] == 'gene_symbol' else True)
        mapping['mapped_id'] = mapping[meta['id_type']]
        clusters = sorted(data.obs.cluster.astype(str).unique())
    finally:
        data.file.close()
    background = sorted(mapping.loc[mapping.included, 'mapped_id'])
    if not background:
        raise InputError('ENRICHMENT_BACKGROUND: no unambiguous detected genes')
    background_set = set(background)
    sets = {term: genes & background_set for term, genes in sets.items()}
    sets = {term: genes for term, genes in sets.items() if min_size <= len(genes) <= max_size}
    if not sets:
        raise InputError('ENRICHMENT_SETS: no gene sets pass background overlap and size filters')
    # Inert lowercase tokens prevent GSEApy's implicit symbol uppercasing.
    tokens = {gene: f'gene_{i:08d}' for i, gene in enumerate(background)}
    token_sets = {term: [tokens[gene] for gene in sorted(genes)] for term, genes in sets.items()}
    markers = pd.read_csv(run / 'markers.csv', dtype={'cluster_id': str, 'gene_id': str})
    mapped = mapping.set_index('gene_id')
    selections, rows = [], []
    for cluster in clusters:
        selected = markers.loc[(markers.cluster_id == cluster) & (markers.logfoldchange > 0) & (markers.pval_adj <= marker_fdr)].copy()
        selected['mapped_id'] = selected.gene_id.map(mapped.mapped_id)
        selected['included'] = selected.gene_id.map(mapped.included).fillna(False).astype(bool)
        selections.append(selected)
        query = sorted(set(selected.loc[selected.included, 'mapped_id']) & background_set)
        if query:
            result = gp.enrich(gene_list=[tokens[g] for g in query], gene_sets=token_sets,
                               background=list(tokens.values()), outdir=None, no_plot=True, verbose=False)
            hits = result.results if isinstance(result.results, pd.DataFrame) else pd.DataFrame()
        else:
            hits = pd.DataFrame()
        by_term = hits.set_index('Term') if not hits.empty else pd.DataFrame()
        for term, genes in sets.items():
            hit = by_term.loc[term] if term in by_term.index else None
            overlap = sorted(genes & set(query))
            rows.append({'cluster_id': cluster, 'term': term, 'query_genes': len(query), 'background_genes': len(background),
                'set_genes_in_background': len(genes), 'overlap_genes': len(overlap),
                'pvalue': float(hit['P-value']) if hit is not None else 1.0, 'genes': ';'.join(overlap)})
    results = pd.DataFrame(rows)
    results['pval_adj'] = false_discovery_control(results.pvalue.to_numpy(), method='bh')
    results = results.sort_values(['pval_adj', 'pvalue', 'cluster_id', 'term'])
    if any(sha256(path) != inputs[key] for key, path in {'processed': run / 'processed.h5ad', 'markers': run / 'markers.csv',
           'manifest': run / 'run_manifest.json', 'gmt': gmt, 'metadata': metadata}.items()):
        raise InputError('ENRICHMENT_STALE: inputs changed during calculation')
    out.mkdir(parents=True)
    results.to_csv(out / 'enrichment.csv', index=False)
    mapping.to_csv(out / 'gene_mapping.csv', index=False)
    pd.concat(selections, ignore_index=True).to_csv(out / 'selected_markers.csv', index=False)
    pd.Series(background, name='gene').to_csv(out / 'background_genes.csv', index=False)
    (out / 'effective_gene_sets.gmt').write_text(''.join(term + '\tbackground-filtered\t' + '\t'.join(sorted(genes)) + '\n' for term, genes in sets.items()), encoding='utf-8')
    summary = {'status': 'completed', 'source_run': str(run), 'input_sha256': inputs, 'gene_sets': meta,
        'gseapy_version': importlib.metadata.version('gseapy'), 'method': 'GSEApy offline over-representation; positive markers from saved top-30 table',
        'marker_fdr': marker_fdr, 'min_set_size': min_size, 'max_set_size': max_size, 'background_size': len(background),
        'tested_sets': len(sets), 'cluster_term_tests': len(results), 'fdr_scope': 'BH across ALL cluster-term tests, including zero-hit p=1 terms',
        'warnings': ['Exploratory enrichment of truncated cluster-vs-rest marker lists; not embryo-replicate differential expression.',
                     'Provided symbols are matched exactly; ambiguous symbols excluded. Genome/annotation version may be unknown.'],
        'significant_tests_fdr_005': int((results.pval_adj <= .05).sum())}
    summary['files'] = ['enrichment.csv', 'gene_mapping.csv', 'selected_markers.csv', 'background_genes.csv',
                        'effective_gene_sets.gmt', 'enrichment_manifest.json', 'report.html']
    (out / 'enrichment_manifest.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')
    html = "<!doctype html><meta charset='utf-8'><title>Exploratory enrichment / 探索性富集</title><h1>Exploratory enrichment / 探索性富集</h1>"
    html += '<p>候选功能解释，不代表已确认细胞身份或胚胎级差异。Exploratory associations; no confirmed cell identity or embryo-level inference.</p>'
    html += '<pre>' + escape(json.dumps(summary, indent=2, ensure_ascii=False)) + '</pre>' + results.head(100).to_html(index=False, escape=True)
    html += '<p>' + ' '.join(f"<a href='{name}'>{name}</a>" for name in summary['files'] if name != 'report.html') + '</p>'
    (out / 'report.html').write_text(html, encoding='utf-8')
    with zipfile.ZipFile(out / 'enrichment_bundle.zip', 'w', compression=zipfile.ZIP_DEFLATED) as bundle:
        for name in summary['files']:
            bundle.write(out / name, arcname=name)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(prog='embryoomics-enrich')
    for key in ('run', 'gmt', 'metadata', 'out'):
        parser.add_argument('--' + key, required=True)
    parser.add_argument('--marker-fdr', type=float, default=.05)
    parser.add_argument('--min-size', type=int, default=5)
    parser.add_argument('--max-size', type=int, default=500)
    args = parser.parse_args(argv)
    try:
        print(json.dumps(enrich_markers(args.run, args.gmt, args.metadata, args.out,
              marker_fdr=args.marker_fdr, min_size=args.min_size, max_size=args.max_size), indent=2))
        return 0
    except (InputError, OSError, KeyError, ValueError) as exc:
        print(str(exc))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
