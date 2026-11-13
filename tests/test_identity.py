import pandas as pd
import pytest
from embryoomics.identity import validate_mapping
from embryoomics.config import InputError


def samples():
    return pd.DataFrame([
        dict(record_id='rna', dataset_id='A', sample_id='s1', species='mus_musculus', stage='E7.0', modality='rna', assay='rna_10x', measurement_unit='cell', cell_id='c1', embryo_id='e1'),
        dict(record_id='atac', dataset_id='A', sample_id='s1', species='mus_musculus', stage='E7.0', modality='atac', assay='atac_fragments', measurement_unit='cell', cell_id='c1', embryo_id='e1')])


def pairs(level='same_cell'):
    return pd.DataFrame([dict(left_record='rna', right_record='atac', pairing_level=level, evidence_url='https://example.invalid/synthetic')])


@pytest.mark.parametrize('level', ['same_cell', 'same_embryo', 'same_sample', 'stage_only', 'unknown', 'unpaired'])
def test_declared_levels_and_audit_outputs(tmp_path, level):
    summary = validate_mapping(samples(), pairs(level), tmp_path / 'mapping')
    assert summary['declared_pairs'] == 1
    assert (tmp_path / 'mapping/sample_mapping.csv').is_file()
    assert summary['unmatched_records'] == (2 if level in ('unknown', 'unpaired') else 0)


@pytest.mark.parametrize('field,value,error', [('dataset_id', 'B', 'PAIR_NAMESPACE'), ('cell_id', '', 'PAIR_ID'),
    ('species', 'homo_sapiens', 'PAIR_CONFLICT'), ('stage', 'E7.5', 'PAIR_CONFLICT')])
def test_false_pairing_rejected(tmp_path, field, value, error):
    data = samples()
    data.loc[1, field] = value
    with pytest.raises(InputError, match=error):
        validate_mapping(data, pairs(), tmp_path / 'invalid')


def test_stage_only_not_upgraded_and_missing_evidence_rejected(tmp_path):
    data = samples()
    data.loc[1, 'dataset_id'] = 'B'
    result = validate_mapping(data, pairs('stage_only'), tmp_path / 'stage')
    assert result['levels'] == {'stage_only': 1}
    links = pairs()
    links.loc[0, 'evidence_url'] = ''
    with pytest.raises(InputError, match='PAIR_EVIDENCE'):
        validate_mapping(samples(), links, tmp_path / 'no-evidence')
