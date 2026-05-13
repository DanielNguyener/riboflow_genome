#!/usr/bin/env python3
"""Replacement for ``rfc stats-percentage`` that operates on the new pipeline's
row-per-stat CSV schema and inserts reads-based retention percentages.

Input CSV layout (from ``rfc merge overall-stats``):

    ,sample1.1,sample1.2,sample2.1,...
    total_reads,...
    clipped_reads,...
    filtered_out,...
    filter_kept,...
    genome_aligned_once,...
    genome_aligned_many,...
    genome_unaligned,...
    genome_primary_alignments,...
    genome_secondary_alignments,...
    genome_total_alignments,...
    qpass_primary_alignments,...
    qpass_secondary_alignments,...
    qpass_total_alignments,...
    dedup_primary_alignments,...
    dedup_secondary_alignments,...
    dedup_total_alignments,...

Output CSV preserves all raw rows and inserts derived rows immediately after
their numerator row:

- ``<numerator>_%`` for between-step retention pairs.
- ``<step>_pct_primary`` for within-step primary-of-total proportions.

The percentages are computed once at the merged level (summed across lanes
where applicable). Computing them per-lane and then summing would be wrong,
which is why the per-lane CSV emitted by ``individual_genome_alignment_stats``
intentionally contains only raw counts.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np


# Each entry is (numerator_row, denominator_row, derived_row_name).
# Order matters: rows are emitted in this sequence, and any row referenced as
# a numerator OR denominator must already exist in the input CSV.
PERCENTAGE_ROWS = [
    ("clipped_reads",              "total_reads",                "clipped_reads_%"),
    ("filtered_out",               "clipped_reads",              "filtered_out_%"),
    ("filter_kept",                "clipped_reads",              "filter_kept_%"),
    ("genome_primary_alignments",  "filter_kept",                "genome_primary_alignments_%"),
    ("genome_primary_alignments",  "genome_total_alignments",    "genome_pct_primary"),
    ("qpass_primary_alignments",   "genome_primary_alignments",  "qpass_primary_alignments_%"),
    ("qpass_primary_alignments",   "qpass_total_alignments",     "qpass_pct_primary"),
    ("dedup_primary_alignments",   "qpass_primary_alignments",   "dedup_primary_alignments_%"),
    ("dedup_primary_alignments",   "dedup_total_alignments",     "dedup_pct_primary"),
]


# Canonical row order in the output. Derived rows are slotted immediately
# after the row that drives them, except for the within-step pct_primary
# rows which are placed after each step's total.
OUTPUT_ROW_ORDER = [
    "total_reads",
    "clipped_reads",
    "clipped_reads_%",
    "filtered_out",
    "filtered_out_%",
    "filter_kept",
    "filter_kept_%",
    "genome_aligned_once",
    "genome_aligned_many",
    "genome_unaligned",
    "genome_primary_alignments",
    "genome_primary_alignments_%",
    "genome_secondary_alignments",
    "genome_total_alignments",
    "genome_pct_primary",
    "qpass_primary_alignments",
    "qpass_primary_alignments_%",
    "qpass_secondary_alignments",
    "qpass_total_alignments",
    "qpass_pct_primary",
    "dedup_primary_alignments",
    "dedup_primary_alignments_%",
    "dedup_secondary_alignments",
    "dedup_total_alignments",
    "dedup_pct_primary",
]


def _safe_pct(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    """Return ``100 * num / denom`` rounded to 2 dp; 0 where denom == 0.

    Division by zero is treated as 0% (not NaN/inf) so downstream consumers
    don't have to special-case empty samples or the no-input case.
    """
    num = pd.to_numeric(numerator, errors="coerce").fillna(0)
    den = pd.to_numeric(denominator, errors="coerce").fillna(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(den == 0, 0.0, 100.0 * num / den)
    return pd.Series(np.round(pct, 2), index=numerator.index)


def stats_percentage(input_csv: Path, output_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(input_csv, header=0, index_col=0)

    for num_row, denom_row, derived_row in PERCENTAGE_ROWS:
        missing = [r for r in (num_row, denom_row) if r not in df.index]
        if missing:
            raise KeyError(
                f"Cannot compute {derived_row}: missing row(s) {missing} in "
                f"{input_csv}. Found rows: {list(df.index)}"
            )
        df.loc[derived_row] = _safe_pct(df.loc[num_row], df.loc[denom_row])

    # Reindex to the canonical output order, keeping any extra rows (e.g.
    # genome_aligned_once) at the bottom so we don't silently drop data.
    extras = [r for r in df.index if r not in OUTPUT_ROW_ORDER]
    final_order = [r for r in OUTPUT_ROW_ORDER if r in df.index] + extras
    df = df.loc[final_order]
    df.to_csv(output_csv)
    return df


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="Input merged-overall-stats CSV (row-per-stat).")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Output CSV with inserted percentage rows.")
    args = parser.parse_args(argv)

    stats_percentage(args.input, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
