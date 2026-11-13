"""Command line entrypoint."""

import argparse
import json
import sys

from .config import InputError
from .pipeline import inspect_input, inspect_doublets, run_analysis
from pathlib import Path
import pandas as pd


def main(argv=None):
    parser = argparse.ArgumentParser(prog="embryoomics")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "run", "preview-doublets"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--input", required=True, help="10x folder, .h5ad, or declared gene-by-cell TSV")
        cmd.add_argument("--config", required=True)
        cmd.add_argument("--metadata", help="Optional metadata.csv with exact cell IDs")
        if name == 'preview-doublets':
            cmd.add_argument('--out', required=True, help='New directory for scores and preview summary')
        if name == "run":
            cmd.add_argument("--out", required=True)
            cmd.add_argument("--confirm-qc", action="store_true", help="Confirm reviewed QC thresholds in config")
            cmd.add_argument('--confirm-doublets', action='store_true', help='Explicitly confirm filtering the reviewed Scrublet candidates')
            cmd.add_argument("--annotations", help="Human-reviewed annotation CSV")
            cmd.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            info = inspect_input(args.input, args.config, args.metadata)
            print(json.dumps({key: info[key] for key in ("cells", "genes", "mt", "warnings", "retained", "files")}, ensure_ascii=False, indent=2))
            print(info["preview"].to_string(index=False))
        elif args.command == 'preview-doublets':
            out = Path(args.out)
            if out.exists():
                raise InputError('OUTPUT_EXISTS: choose a new preview directory')
            preview = inspect_doublets(args.input, args.config, args.metadata)
            out.mkdir(parents=True)
            preview['audit'].to_csv(out / 'doublet_scores.csv')
            pd.DataFrame({'doublet_score': preview['simulated']}).to_csv(out / 'doublet_simulated_scores.csv', index=False)
            (out / 'doublet_preview.json').write_text(json.dumps({'summary': preview['summary'], 'retained': preview['retained'], 'input_identity': preview['input_identity'], 'config': preview['config']}, indent=2), encoding='utf-8')
            print(json.dumps(preview['summary'], indent=2))
        else:
            if not args.confirm_qc:
                raise InputError("QC_CONFIRMATION: run validate, review its QC preview, then pass --confirm-qc")
            print(json.dumps(run_analysis(args.config, args.input, args.out, args.metadata, args.annotations, args.overwrite,
                                          confirm_doublets=args.confirm_doublets), ensure_ascii=False, indent=2))
        return 0
    except (InputError, FileNotFoundError) as exc:
        print(str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
