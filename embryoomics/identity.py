"""Validate declared sample identities and evidence-backed pairing levels."""
import argparse
import json
from pathlib import Path
import pandas as pd
from .config import InputError


ASSAYS = {'rna_10x': ('rna', 'cell'), 'mass_spec_quantification': ('protein', 'sample'),
          'atac_fragments': ('atac', 'cell')}
LEVELS = ('same_cell', 'same_embryo', 'same_sample', 'stage_only', 'unpaired', 'unknown')


def validate_mapping(samples, pairs, out):
    required = {'record_id', 'dataset_id', 'sample_id', 'species', 'stage', 'modality', 'assay', 'measurement_unit'}
    if not required.issubset(samples) or samples[list(required)].eq('').any().any():
        raise InputError('IDENTITY_FIELDS: missing or empty required sample fields')
    if samples.record_id.duplicated().any() or samples.duplicated(['dataset_id', 'sample_id', 'modality', 'cell_id'] if 'cell_id' in samples else ['dataset_id', 'sample_id', 'modality']).any():
        raise InputError('IDENTITY_DUPLICATE: record and composite sample/cell keys must be unique')
    for row in samples.itertuples():
        if row.species not in ('mus_musculus', 'homo_sapiens') or row.assay not in ASSAYS or ASSAYS[row.assay] != (row.modality, row.measurement_unit):
            raise InputError('IDENTITY_ASSAY: unsupported species or inconsistent assay/modality/measurement unit')
    if not {'left_record', 'right_record', 'pairing_level', 'evidence_url'}.issubset(pairs):
        raise InputError('PAIR_FIELDS: declare both endpoints, pairing_level and evidence_url')
    if pairs.duplicated(['left_record', 'right_record']).any():
        raise InputError('PAIR_DUPLICATE: duplicate endpoint pair')
    lookup = samples.set_index('record_id')
    audited = []
    for row in pairs.itertuples():
        if row.left_record == row.right_record or row.left_record not in lookup.index or row.right_record not in lookup.index:
            raise InputError('PAIR_ENDPOINT: missing or identical endpoints')
        left, right = lookup.loc[row.left_record], lookup.loc[row.right_record]
        if row.pairing_level not in LEVELS:
            raise InputError('PAIR_LEVEL: unsupported pairing level')
        level = row.pairing_level
        if level not in ('unpaired', 'unknown'):
            if not str(row.evidence_url).startswith(('https://', 'http://', 'doi:')):
                raise InputError('PAIR_EVIDENCE: an explicit source reference is required')
            if left.species != right.species or left.stage != right.stage or left.stage.lower() in ('unknown', 'unassigned'):
                raise InputError('PAIR_CONFLICT: species/stage mismatch or unknown stage')
            if left.modality == right.modality:
                raise InputError('PAIR_MODALITY: cross-modal pair must have different modalities')
        if level in ('same_cell', 'same_embryo', 'same_sample'):
            if left.dataset_id != right.dataset_id:
                raise InputError('PAIR_NAMESPACE: same-named IDs from different datasets cannot imply pairing')
            for optional in ('embryo_id', 'replicate_id'):
                if optional in samples and left[optional] not in ('', 'unknown') and right[optional] not in ('', 'unknown') and left[optional] != right[optional]:
                    raise InputError('PAIR_CONFLICT: declared embryo or replicate IDs disagree')
            key = {'same_cell': 'cell_id', 'same_embryo': 'embryo_id', 'same_sample': 'sample_id'}[level]
            if key not in samples or not left[key] or left[key].lower() == 'unknown' or left[key] != right[key]:
                raise InputError('PAIR_ID: matching documented IDs are required for the declared level')
            if level == 'same_cell' and (left.measurement_unit != 'cell' or right.measurement_unit != 'cell' or left.sample_id != right.sample_id):
                raise InputError('PAIR_UNIT: same-cell pairing requires cell measurements from the same sample')
        audited.append({'left_record': row.left_record, 'right_record': row.right_record,
                        'pairing_level': level, 'evidence_url': row.evidence_url, 'status': 'declared_not_independently_verified'})
    out = Path(out)
    if out.exists():
        raise InputError('OUTPUT_EXISTS: choose a new mapping directory')
    out.mkdir(parents=True)
    pd.DataFrame(audited, columns=['left_record', 'right_record', 'pairing_level', 'evidence_url', 'status']).to_csv(out / 'sample_mapping.csv', index=False)
    paired = set(pairs.loc[~pairs.pairing_level.isin(['unknown', 'unpaired']), ['left_record', 'right_record']].to_numpy().ravel())
    unmatched = samples.loc[~samples.record_id.isin(paired)]
    unmatched.to_csv(out / 'unmatched_samples.csv', index=False)
    summary = {'records': len(samples), 'declared_pairs': len(pairs), 'unmatched_records': len(unmatched),
               'levels': pairs.pairing_level.value_counts().to_dict(),
               'warning': 'Validates declared metadata only; does not establish biological pairing or approve integration'}
    (out / 'pairing_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser()
    for key in ('samples', 'pairs', 'out'):
        parser.add_argument('--' + key, required=True)
    args = parser.parse_args(argv)
    try:
        result = validate_mapping(pd.read_csv(args.samples, dtype=str, keep_default_na=False),
                                  pd.read_csv(args.pairs, dtype=str, keep_default_na=False), args.out)
        print(json.dumps(result, indent=2))
        return 0
    except (InputError, OSError, ValueError) as exc:
        print(str(exc))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
