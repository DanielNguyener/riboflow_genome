#!/usr/bin/env python3
"""Parse STAR Log.final.out into a TSV row for downstream aggregation.

Two modes:

  parse_star_log.py <Log.final.out>
      Read a STAR Log.final.out and write a single-row TSV (with header)
      to stdout. Columns are stable so consumers can DictReader the result.

  parse_star_log.py <Log.final.out> --secondary-from-bam <bam>
      Same as above, plus an additional `secondary_alignments` column
      sourced from `samtools view -c -f 256 <bam>`. The count is cached
      as a sidecar file at <bam>.secondary.count and reused on subsequent
      runs (skipped if the sidecar is present and newer than the BAM).

Exit non-zero on any failure to parse a required STAR field — silent
defaults would mask format drift.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

# STAR Log.final.out field labels we care about. STAR uses "|" as the
# key/value separator inside this file.
STAR_FIELD_MAP = {
    "Number of input reads": "total_input_reads",
    "Uniquely mapped reads number": "uniquely_mapped",
    "Uniquely mapped reads %": "uniquely_mapped_pct",
    "Number of reads mapped to multiple loci": "multi_loci_mapped",
    "% of reads mapped to multiple loci": "multi_loci_mapped_pct",
    "Number of reads mapped to too many loci": "too_many_loci",
    "% of reads mapped to too many loci": "too_many_loci_pct",
    "Number of reads unmapped: too many mismatches": "unmapped_too_many_mismatches",
    "Number of reads unmapped: too short": "unmapped_too_short",
    "Number of reads unmapped: other": "unmapped_other",
    "Mismatch rate per base, %": "mismatch_rate_per_base_pct",
    "Mapping speed, Million of reads per hour": "mapping_speed_million_per_hour",
}

OUTPUT_COLUMNS = [
    "total_input_reads",
    "uniquely_mapped",
    "uniquely_mapped_pct",
    "multi_loci_mapped",
    "multi_loci_mapped_pct",
    "too_many_loci",
    "too_many_loci_pct",
    "unmapped_too_many_mismatches",
    "unmapped_too_short",
    "unmapped_other",
    "unmapped_total",
    "primary_aligned_total",
    "mismatch_rate_per_base_pct",
    "mapping_speed_million_per_hour",
]


def parse_star_log(path: Path) -> dict[str, str]:
    """Return a dict keyed by OUTPUT_COLUMNS. Values stay as strings so
    the caller can decide whether to coerce; STAR's percentage fields
    carry a trailing '%' which we strip but keep numeric-as-string."""
    fields: dict[str, str] = {}
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if "|" not in line:
                continue
            label, value = line.split("|", 1)
            label = label.strip()
            value = value.strip().rstrip("%").strip()
            if label in STAR_FIELD_MAP:
                fields[STAR_FIELD_MAP[label]] = value

    # Required fields STAR always emits — bail loudly if missing.
    required = ("total_input_reads", "uniquely_mapped", "multi_loci_mapped")
    missing = [r for r in required if r not in fields]
    if missing:
        raise ValueError(
            f"STAR log {path} is missing required fields: {', '.join(missing)}"
        )

    # Derived fields.
    unmapped_total = (
        int(fields.get("unmapped_too_many_mismatches", 0) or 0)
        + int(fields.get("unmapped_too_short", 0) or 0)
        + int(fields.get("unmapped_other", 0) or 0)
    )
    primary_total = int(fields["uniquely_mapped"]) + int(fields["multi_loci_mapped"])

    fields["unmapped_total"] = str(unmapped_total)
    fields["primary_aligned_total"] = str(primary_total)

    # Fill any optional column that STAR didn't emit so the TSV has a
    # consistent shape across versions.
    for col in OUTPUT_COLUMNS:
        fields.setdefault(col, "")

    return fields


def secondary_count(bam_path: Path) -> int:
    """Return the count of secondary alignments (FLAG 256) in `bam_path`,
    using a `<bam>.secondary.count` sidecar file as a cache."""
    sidecar = bam_path.with_suffix(bam_path.suffix + ".secondary.count")
    if sidecar.exists() and sidecar.stat().st_mtime >= bam_path.stat().st_mtime:
        return int(sidecar.read_text().strip())

    # samtools view -c -f 256 prints just the count to stdout.
    result = subprocess.run(
        ["samtools", "view", "-c", "-f", "256", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    count = int(result.stdout.strip())
    sidecar.write_text(f"{count}\n")
    return count


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log_final_out", type=Path, help="STAR Log.final.out path")
    parser.add_argument(
        "--secondary-from-bam",
        type=Path,
        default=None,
        help=(
            "Optional BAM path. If given, append a `secondary_alignments` "
            "column counted via `samtools view -c -f 256`, cached as "
            "<bam>.secondary.count."
        ),
    )
    args = parser.parse_args(argv)

    if not args.log_final_out.is_file():
        print(f"ERROR: STAR log not found: {args.log_final_out}", file=sys.stderr)
        return 2

    fields = parse_star_log(args.log_final_out)

    columns = list(OUTPUT_COLUMNS)
    if args.secondary_from_bam is not None:
        if not args.secondary_from_bam.is_file():
            print(
                f"ERROR: BAM not found: {args.secondary_from_bam}",
                file=sys.stderr,
            )
            return 2
        fields["secondary_alignments"] = str(secondary_count(args.secondary_from_bam))
        columns.append("secondary_alignments")

    sys.stdout.write("\t".join(columns) + "\n")
    sys.stdout.write("\t".join(fields[c] for c in columns) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
