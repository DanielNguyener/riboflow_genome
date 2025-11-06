#!/usr/bin/env python3

"""
Update RNA-seq merged statistics CSV with dedup count.
This script reads a temporary RNA-seq merged stats CSV and updates it with
actual merged dedup count from deduplication.
"""

import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(
        description='Update RNA-seq merged stats CSV with dedup count'
    )
    parser.add_argument(
        '--dedup-count-file',
        required=True,
        help='Path to merged dedup count file'
    )
    parser.add_argument(
        '--input-csv',
        required=True,
        help='Path to input temporary CSV file'
    )
    parser.add_argument(
        '--output-csv',
        required=True,
        help='Path to output CSV file'
    )
    
    args = parser.parse_args()
    
    # Read the merged dedup count
    try:
        with open(args.dedup_count_file, 'r') as f:
            merged_dedup_count = f.read().strip()
    except FileNotFoundError:
        print(f"Error: Dedup count file not found: {args.dedup_count_file}", file=sys.stderr)
        sys.exit(1)
    
    # Read the merged stats CSV
    try:
        with open(args.input_csv, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = reader.fieldnames
    except FileNotFoundError:
        print(f"Error: Input CSV file not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)
    
    # Update genome_after_dedup with merged count
    for row in rows:
        if 'genome_after_dedup' in row:
            row['genome_after_dedup'] = merged_dedup_count
    
    # Write updated CSV
    with open(args.output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()

