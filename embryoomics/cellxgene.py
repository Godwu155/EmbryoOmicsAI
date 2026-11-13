"""Optional CELLxGENE exchange without importing its incompatible runtime."""

import argparse
import hashlib
import json
from io import StringIO
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .config import InputError
from .evidence import validate_annotations
from .review import save_annotation_review


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def export_view(run, out):
    run, out = Path(run).resolve(), Path(out).resolve()
    if out.exists():
        raise InputError('OUTPUT_EXISTS: choose a new viewer export directory')
    source = run / 'processed.h5ad'
    annotations_path = run / 'annotations.csv'
    before = {'processed_sha256': sha256(source), 'annotations_sha256': sha256(annotations_path)}
    data = ad.read_h5ad(source, backed='r')
    try:
        if 'X_umap' not in data.obsm or 'cluster' not in data.obs:
            raise InputError('VIEWER_SCHEMA: processed data needs cluster and X_umap')
        annotations = pd.read_csv(annotations_path, dtype=str, keep_default_na=False)
        annotations = validate_annotations(annotations, data.obs.cluster.astype(str).unique())
        obs = data.obs[[key for key in ('cluster', 'stage', 'sample_id', 'embryo_id', 'batch', 'doublet_score', 'predicted_doublet') if key in data.obs]].copy()
        for key in ('cluster', 'stage', 'sample_id', 'embryo_id', 'batch'):
            if key in obs:
                obs[key] = obs[key].astype(str).astype('category')
        var = pd.DataFrame(index=data.var_names)
        var['gene_id'] = data.var_names.astype(str)
        var['gene_symbol'] = data.var.get('gene_symbol', var.gene_id).astype('string').fillna('').astype(str)
        var['viewer_gene'] = [f'{symbol or gene} [{gene}]' for gene, symbol in zip(var.gene_id, var.gene_symbol)]
        if not var.viewer_gene.is_unique:
            raise InputError('VIEWER_GENES: display labels must be unique')
        matrix = data.X.to_memory() if hasattr(data.X, 'to_memory') else data.X
        view = ad.AnnData(matrix.astype(np.float32).tocsc(), obs=obs, var=var)
        view.obsm['X_umap'] = np.asarray(data.obsm['X_umap']).copy()
        labels = pd.DataFrame({'candidate_cell_type': obs.cluster.astype(str).map(annotations.set_index('cluster_id').candidate_cell_type)}, index=obs.index)
        labels.index.name = 'cell_id'
        if sha256(source) != before['processed_sha256'] or sha256(annotations_path) != before['annotations_sha256']:
            raise InputError('VIEWER_STALE: source changed during export; export again')
        out.mkdir(parents=True)
        view.write_h5ad(out / 'viewer.h5ad', compression='gzip')
        labels.to_csv(out / 'cell_annotations.csv')
        var.to_csv(out / 'gene_map.csv', index=False)
        provenance = {'source_run': str(run), **before, 'cells': view.n_obs, 'genes': view.n_vars,
            'viewer_sha256': sha256(out / 'viewer.h5ad'), 'expression': 'X log1p normalized expression, float32 CSC browsing copy',
            'embedding': 'Existing X_umap, no recomputation', 'cellxgene_version': '1.3.0',
            'annotation_column': 'candidate_cell_type', 'notes': 'User cell labels remain candidates; source run is unchanged'}
        (out / 'exchange.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
        return provenance
    finally:
        data.file.close()


def import_proposals(exchange, labels_path, out):
    """Propose only unanimous per-cluster labels; preserve evidence and flag splits."""
    exchange, out = Path(exchange).resolve(), Path(out).resolve()
    provenance = json.loads((exchange / 'exchange.json').read_text(encoding='utf-8'))
    run = Path(provenance['source_run'])
    if out.exists():
        raise InputError('OUTPUT_EXISTS: choose a new proposal directory')
    if sha256(run / 'processed.h5ad') != provenance['processed_sha256'] or sha256(run / 'annotations.csv') != provenance['annotations_sha256']:
        raise InputError('VIEWER_STALE: source run or annotation review changed; export a fresh view')
    if sha256(exchange / 'viewer.h5ad') != provenance['viewer_sha256']:
        raise InputError('VIEWER_STALE: viewer data changed since export')
    labels_version = sha256(labels_path)
    content = Path(labels_path).read_text(encoding='utf-8-sig')
    while content.startswith('#'):
        _, _, content = content.partition('\n')
    labels = pd.read_csv(StringIO(content), dtype=str, keep_default_na=False)
    if len(labels.columns) < 2:
        raise InputError('VIEWER_LABELS: expected a cell ID index and candidate_cell_type column')
    labels = labels.rename(columns={labels.columns[0]: 'cell_id'})
    if 'candidate_cell_type' not in labels or labels.cell_id.duplicated().any():
        raise InputError('VIEWER_LABELS: unique cell IDs and candidate_cell_type are required')
    data = ad.read_h5ad(run / 'processed.h5ad', backed='r')
    try:
        if set(labels.cell_id) != set(data.obs_names):
            raise InputError('VIEWER_LABELS: cell IDs must exactly match the exported run')
        obs = data.obs[['cluster']].copy()
    finally:
        data.file.close()
    labels = labels.set_index('cell_id').loc[obs.index, 'candidate_cell_type'].str.strip()
    proposals = pd.read_csv(run / 'annotations.csv', dtype=str, keep_default_na=False)
    proposals = validate_annotations(proposals, obs.cluster.astype(str).unique())
    proposals['review_status'] = 'in_review'
    conflicts = []
    for cluster, cells in obs.groupby('cluster', observed=True):
        values = labels.loc[cells.index]
        counts = values.value_counts().to_dict()
        unanimous = len(counts) == 1 and values.iloc[0] not in ('', 'unassigned', 'unknown')
        row = proposals.cluster_id == str(cluster)
        if unanimous:
            proposals.loc[row, 'candidate_cell_type'] = values.iloc[0]
        else:
            conflicts.append({'cluster_id': str(cluster), 'cells': len(cells), 'labels': json.dumps(counts, ensure_ascii=False)})
        prefix = proposals.loc[row, 'reviewer_note'].iloc[0]
        note = f'CELLxGENE cell-label proposal: {len(cells)} cells; ' + ('unanimous candidate, verify evidence' if unanimous else 'unresolved or mixed labels; previous candidate retained')
        proposals.loc[row, 'reviewer_note'] = (prefix + '\n' + note).strip()
    # Recheck the source after reading before writing review artifacts.
    if sha256(run / 'annotations.csv') != provenance['annotations_sha256'] or sha256(run / 'processed.h5ad') != provenance['processed_sha256']:
        raise InputError('VIEWER_STALE: annotations changed while preparing proposals')
    if sha256(labels_path) != labels_version:
        raise InputError('VIEWER_STALE: CELLxGENE labels changed while preparing proposals; wait for All saved and retry')
    out.mkdir(parents=True)
    proposals.to_csv(out / 'annotation_proposals.csv', index=False)
    pd.DataFrame(conflicts, columns=['cluster_id', 'cells', 'labels']).to_csv(out / 'annotation_conflicts.csv', index=False)
    summary = {'clusters': len(proposals), 'unresolved_clusters': len(conflicts),
        'labels_sha256': labels_version, 'source_run': str(run),
        'source_annotations_sha256': provenance['annotations_sha256'], 'source_modified': False}
    summary['source_processed_sha256'] = provenance['processed_sha256']
    (out / 'proposal_manifest.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def apply_proposals(proposal_dir, *, confirm_reviewed=False):
    if not confirm_reviewed:
        raise InputError('VIEWER_CONFIRMATION: review proposal and conflict CSVs, then pass --confirm-reviewed')
    proposal_dir = Path(proposal_dir)
    manifest = json.loads((proposal_dir / 'proposal_manifest.json').read_text(encoding='utf-8'))
    run = Path(manifest['source_run'])
    if sha256(run / 'processed.h5ad') != manifest['source_processed_sha256']:
        raise InputError('VIEWER_STALE: analysis changed after proposal generation')
    proposals = pd.read_csv(proposal_dir / 'annotation_proposals.csv', dtype=str, keep_default_na=False)
    # Reuse the normal review lock, stale-version guard, history and report/ZIP refresh.
    return save_annotation_review(run, proposals, manifest['source_annotations_sha256'])


def main(argv=None):
    parser = argparse.ArgumentParser(prog='embryoomics-cellxgene')
    sub = parser.add_subparsers(dest='command', required=True)
    export = sub.add_parser('export')
    export.add_argument('--run', required=True)
    export.add_argument('--out', required=True)
    propose = sub.add_parser('propose')
    propose.add_argument('--exchange', required=True)
    propose.add_argument('--labels', required=True)
    propose.add_argument('--out', required=True)
    apply = sub.add_parser('apply')
    apply.add_argument('--proposals', required=True)
    apply.add_argument('--confirm-reviewed', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.command == 'export':
            result = export_view(args.run, args.out)
        elif args.command == 'propose':
            result = import_proposals(args.exchange, args.labels, args.out)
        else:
            result = apply_proposals(args.proposals, confirm_reviewed=args.confirm_reviewed)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except (InputError, OSError, KeyError, ValueError) as exc:
        print(str(exc))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
