import json
import pytest
from embryoomics.pipeline import run_analysis
from test_contract import sample


def test_output_failure_records_failed_status_without_overwriting_other_runs(sample, tmp_path, monkeypatch):
    import embryoomics.pipeline as pipeline
    folder, config_path, _, _ = sample
    out = tmp_path / 'failed'
    def fail(*args, **kwargs):
        raise OSError('test output failure')
    monkeypatch.setattr(pipeline, 'write_report', fail)
    with pytest.raises(OSError, match='test output failure'):
        run_analysis(config_path, folder, out)
    assert json.loads((out / 'run_status.json').read_text())['status'] == 'failed'
    assert 'test output failure' in (out / 'run.log').read_text()
    before = (out / 'run_status.json').read_bytes()
    with pytest.raises(ValueError, match='OUTPUT_EXISTS'):
        run_analysis(config_path, folder, out)
    assert (out / 'run_status.json').read_bytes() == before
