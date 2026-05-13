#!/usr/bin/env python3
"""RNA-seq mirror of ``update_merged_stats_with_counts.py``.

Overrides the three ``dedup_*_alignments`` rows in the per-sample RNA-seq
merged stats CSV with counts derived from ``samtools view -c`` on the
sample-level post-dedup BAM. The row-per-stat layout matches the ribo-seq
side; this script is a thin re-export to keep the call sites symmetric.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def _read_first_int(path: Path) -> str:
    with path.open() as fh:
        return fh.read().strip().split()[0]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dedup-total-file",     required=True, type=Path)
    parser.add_argument("--dedup-primary-file",   required=True, type=Path)
    parser.add_argument("--dedup-secondary-file", required=True, type=Path)
    parser.add_argument("--input-csv",  required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    args = parser.parse_args(argv)

    overrides = {
        "dedup_total_alignments":     _read_first_int(args.dedup_total_file),
        "dedup_primary_alignments":   _read_first_int(args.dedup_primary_file),
        "dedup_secondary_alignments": _read_first_int(args.dedup_secondary_file),
    }

    if not args.input_csv.exists():
        print(f"Error: input CSV {args.input_csv} not found", file=sys.stderr)
        return 1

    with args.input_csv.open(newline="") as fh:
        rows = list(csv.reader(fh))

    for row in rows[1:]:
        if row and len(row) >= 2 and row[0] in overrides:
            row[1] = overrides[row[0]]

    with args.output_csv.open("w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
