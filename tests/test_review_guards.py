from pathlib import Path
import pytest
import yaml

from embryoomics.config import InputError, load_config
from embryoomics.io import review_identity


@pytest.mark.parametrize('field,value', [('leiden_resolution', float('nan')), ('leiden_resolution', float('inf')), ('leiden_resolution', True), ('max_pct_mt', True)])
def test_invalid_numerical_config(tmp_path, field, value):
    source = Path(__file__).resolve().parents[1] / 'examples/GSM8559287.yaml'
    config = yaml.safe_load(source.read_text(encoding='utf-8'))
    config['qc' if field == 'max_pct_mt' else 'analysis'][field] = value
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump(config), encoding='utf-8')
    with pytest.raises(InputError, match='CONFIG'):
        load_config(path)


def test_preview_identity_tracks_config_and_input(tmp_path):
    path = tmp_path / 'counts.tsv'
    path.write_text('Gene\ta\ng1\t1\n', encoding='utf-8')
    config = {'qc': {'min_genes': 200}}
    before = review_identity(path, config)
    assert before != review_identity(path, {'qc': {'min_genes': 500}})
    path.write_text('Gene\ta\ngene1\t123\n', encoding='utf-8')
    assert before != review_identity(path, config)
