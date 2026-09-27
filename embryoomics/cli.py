"""Command line entrypoint."""

import argparse
import json
import sys

from .config import InputError
from .pipeline import inspect_input, run_analysis


def main(argv=None):
    parser = argparse.ArgumentParser(prog="embryoomics")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "run"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--input", required=True, help="10x folder or .h5ad")
        cmd.add_argument("--config", required=True)
        cmd.add_argument("--metadata", help="Optional metadata.csv with exact cell IDs")
        if name == "run":
            cmd.add_argument("--out", required=True)
            cmd.add_argument("--confirm-qc", action="store_true", help="Confirm reviewed QC thresholds in config")
            cmd.add_argument("--annotations", help="Human-reviewed annotation CSV")
            cmd.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            info = inspect_input(args.input, args.config, args.metadata)
            print(json.dumps({key: info[key] for key in ("cells", "genes", "mt", "warnings", "retained", "files")}, ensure_ascii=False, indent=2))
            print(info["preview"].to_string(index=False))
        else:
            if not args.confirm_qc:
                raise InputError("QC_CONFIRMATION: run validate, review its QC preview, then pass --confirm-qc")
            print(json.dumps(run_analysis(args.config, args.input, args.out, args.metadata, args.annotations, args.overwrite), ensure_ascii=False, indent=2))
        return 0
    except (InputError, FileNotFoundError) as exc:
        print(str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
