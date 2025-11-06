#!/usr/bin/env python3

"""
Update merged statistics CSV with P-site count only.
This script reads a temporary merged stats CSV and updates it with actual
merged P-site count from P-site correction.
"""

import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(
        description='Update merged stats CSV with P-site count'
    )
    parser.add_argument(
        '--psite-count-file',
        required=True,
        help='Path to merged P-site count file'
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
    
    # Read the merged P-site count
    try:
        with open(args.psite_count_file, 'r') as f:
            merged_psite_count = f.read().strip()
    except FileNotFoundError:
        print(f"Error: P-site count file not found: {args.psite_count_file}", file=sys.stderr)
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
    
    # Update genome_after_psite with merged count
    for row in rows:
        if 'genome_after_psite' in row:
            row['genome_after_psite'] = merged_psite_count
    
    # Write updated CSV
    with open(args.output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()

