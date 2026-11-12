"""Persist annotation edits and regenerate reports without rerunning analysis."""

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import uuid
import zipfile

import pandas as pd

from .config import InputError
from .evidence import validate_annotations
from .reporting import write_report


REVIEW_STATUSES = ['unreviewed', 'in_review', 'reviewed', 'rejected']


def annotation_version(run):
    return hashlib.sha256((Path(run) / 'annotations.csv').read_bytes()).hexdigest()


def save_annotation_review(run, edited, expected_version):
    run = Path(run).resolve()
    lock = run / '.annotation-save.lock'
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise InputError('ANNOTATIONS_BUSY: another save is in progress; try again shortly') from exc
    try:
        os.close(descriptor)
        if annotation_version(run) != expected_version:
            raise InputError('ANNOTATIONS_CONFLICT: annotations changed since this page was loaded; reload before editing')
        if not (run / 'annotation_history').resolve().is_relative_to(run):
            raise InputError('REPORT_PATH: annotation history is outside the selected run')
        manifest = json.loads((run / 'run_manifest.json').read_text(encoding='utf-8'))
        clusters = pd.read_csv(run / 'cluster_composition.csv', dtype=str).cluster.unique().tolist()
        annotations = validate_annotations(edited, clusters)
        if not annotations.review_status.isin(REVIEW_STATUSES).all():
            raise InputError('ANNOTATIONS: invalid review status')
        before = pd.read_csv(run / 'annotations.csv', dtype=str, keep_default_na=False)
        if annotations.equals(validate_annotations(before, clusters)):
            return {'changed': False, 'version': expected_version}
        revision = uuid.uuid4().hex
        timestamp = datetime.now(timezone.utc).isoformat()
        history_name = f'annotation_history/{revision}.json'
        history = {'revision': revision, 'saved_at_utc': timestamp, 'previous_sha256': expected_version,
                   'source': 'local_annotation_editor', 'before': before.to_dict(orient='records'),
                   'after': annotations.to_dict(orient='records')}
        manifest['annotation_review'] = {'revision': revision, 'saved_at_utc': timestamp, 'history_file': history_name}
        manifest['files'] = list(dict.fromkeys([*manifest['files'], history_name]))
        # Stage every changed artifact, including the ZIP, before replacing any live file.
        with tempfile.TemporaryDirectory(prefix='.annotation-save-', dir=run) as temporary:
            stage = Path(temporary)
            (stage / 'annotation_history').mkdir()
            (stage / history_name).write_text(json.dumps(history, ensure_ascii=False, indent=2), encoding='utf-8')
            annotations.to_csv(stage / 'annotations.csv', index=False)
            write_report(stage, manifest, pd.read_csv(run / 'qc_preview.csv'),
                         pd.read_csv(run / 'markers.csv'), annotations)
            with zipfile.ZipFile(stage / 'report_bundle.zip', 'w', compression=zipfile.ZIP_DEFLATED) as bundle:
                for name in manifest['files']:
                    if not isinstance(name, str) or Path(name).is_absolute() or '..' in Path(name).parts:
                        raise InputError('REPORT_PATH: manifest contains an invalid artifact path')
                    source = stage / name if (stage / name).is_file() else run / name
                    if not source.resolve().is_relative_to(run):
                        raise InputError('REPORT_PATH: artifact is outside the selected run')
                    bundle.write(source, arcname=name)
            replacements = [history_name, 'annotations.csv', 'report.html', 'run_manifest.json', 'report_bundle.zip']
            if annotation_version(run) != expected_version:
                raise InputError('ANNOTATIONS_CONFLICT: annotations changed while preparing the report; reload before editing')
            backup = stage / 'backup'
            backup.mkdir()
            replaced = []
            try:
                for name in replacements:
                    target = run / name
                    existed = target.exists()
                    if existed:
                        shutil.copy2(target, backup / Path(name).name)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(stage / name, target)
                    replaced.append((name, existed))
            except OSError:
                for name, existed in reversed(replaced):
                    if existed:
                        os.replace(backup / Path(name).name, run / name)
                    else:
                        (run / name).unlink(missing_ok=True)
                raise
        return {'changed': True, 'version': annotation_version(run), 'revision': revision, 'saved_at_utc': timestamp}
    finally:
        lock.unlink(missing_ok=True)
