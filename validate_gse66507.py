"""Download and inspect the GSE66507 supplementary count matrix.

Run: python validate_gse66507.py --output data/GSE66507_counts.txt.gz
This script reports file-level evidence; it does not perform biological QC.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import urllib.request
from pathlib import Path


URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66507/suppl/"
    "GSE66507_human_blastocyst_rnaseq_counts.txt.gz"
)


def download(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    try:
        with urllib.request.urlopen(URL, timeout=60) as response, temporary.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def inspect(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            digest.update(chunk)

    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            raise ValueError("Expected a gene identifier column and at least one cell column")
        cells = header[1:]
        if any(not cell.strip() for cell in cells) or len(cells) != len(set(cells)):
            raise ValueError("Cell labels are blank or duplicated")

        genes = set()
        rows = 0
        nonzero = 0
        maximum = 0
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Row {line_number}: expected {len(header)} columns, got {len(row)}")
            gene = row[0].strip()
            if not gene or gene in genes:
                raise ValueError(f"Row {line_number}: blank or duplicate gene identifier: {gene!r}")
            genes.add(gene)
            rows += 1
            for value in row[1:]:
                try:
                    count = float(value)
                except ValueError as exc:
                    raise ValueError(f"Row {line_number}: nonnumeric value {value!r}") from exc
                if not math.isfinite(count) or count < 0 or not count.is_integer():
                    raise ValueError(f"Row {line_number}: invalid raw count {value!r}")
                nonzero += count > 0
                maximum = max(maximum, int(count))
        if not rows:
            raise ValueError("No gene rows found")

    return {
        "accession": "GSE66507",
        "source_url": URL,
        "file": str(path.resolve()),
        "compressed_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
        "orientation": "genes_by_cells",
        "genes": rows,
        "cells": len(cells),
        "cell_labels": cells,
        "nonzero_values": nonzero,
        "maximum_count": maximum,
        "nonnegative_integer_values": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/GSE66507_counts.txt.gz"))
    parser.add_argument("--report", type=Path, default=Path("data/GSE66507_validation.json"))
    parser.add_argument("--existing", action="store_true", help="Inspect a previously downloaded file")
    args = parser.parse_args()
    if not args.existing:
        download(args.output)
    report = inspect(args.output)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(report, ensure_ascii=False, indent=2)
    args.report.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
