#!/usr/bin/env python3
"""Override the dedup_* alignment-count rows in a per-sample merged stats CSV.

``rfc sum-stats`` sums per-lane stats into a per-sample CSV. For
``dedup_method`` in {``position``, ``umicollapse``}, the per-lane dedup
counts come from awk-splitting the merged post-dedup BED — which gives the
correct per-lane attribution but loses primary/secondary precision when
summed naively. The authoritative merged counts are the alignment counts
computed by ``samtools view -c -F 2304 | -f 256 | -c`` directly on the
merged post-dedup BAM, which this script substitutes into the row-per-stat
sum-stats output.

The input CSV layout is the row-per-stat format produced by
``rfcommands.sum_stats``::

    ,<sample>
    total_reads,N
    ...
    dedup_total_alignments,N
    dedup_primary_alignments,N
    dedup_secondary_alignments,N
    ...

This script rewrites the three ``dedup_*_alignments`` rows in place and
emits a new CSV otherwise identical to the input.
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

    # First row is the header (`,<sample>`); subsequent rows are
    # (stat_name, value). We rewrite the value column in-place when the stat
    # name matches a key we're overriding.
    for row in rows[1:]:
        if row and len(row) >= 2 and row[0] in overrides:
            row[1] = overrides[row[0]]

    with args.output_csv.open("w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
