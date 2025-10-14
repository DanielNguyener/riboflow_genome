#!/usr/bin/env python3

import argparse
import csv
import sys


def load_offsets(offset_csv, experiment_id):
    offsets = {}
    found_experiment = False

    try:
        with open(offset_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['Experiment'] == experiment_id:
                    found_experiment = True
                    read_length = int(row['Read Length'])
                    offset = int(row['P-site Offset'])
                    offsets[read_length] = offset
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: Offset file not found: {offset_csv}\n")
        sys.exit(1)
    except KeyError as e:
        sys.stderr.write(f"ERROR: Missing column in offset CSV: {e}\n")
        sys.exit(1)

    if not found_experiment:
        sys.stderr.write(f"ERROR: Experiment ID '{experiment_id}' not found in offset file\n")
        sys.exit(1)

    if not offsets:
        sys.stderr.write(f"ERROR: No offsets found for experiment '{experiment_id}'\n")
        sys.exit(1)

    return offsets


def apply_psite_correction_bed(input_bed, output_bed, offsets, sample_name):
    reads_processed = 0
    reads_written = 0
    reads_skipped_no_offset = 0
    reads_skipped_invalid_pos = 0

    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            reads_processed += 1
            fields = line.split('\t')

            if len(fields) < 6:
                sys.stderr.write(f"WARNING: Skipping line with fewer than 6 columns: {line}\n")
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            score = fields[4]
            strand = fields[5]

            read_length = end - start

            if read_length not in offsets:
                reads_skipped_no_offset += 1
                continue

            offset = offsets[read_length]

            if strand == '-':
                psite_end = end - offset
                psite_start = psite_end - 1
            else:
                psite_start = start + offset
                psite_end = psite_start + 1

            if psite_start < 0 or psite_end <= psite_start:
                reads_skipped_invalid_pos += 1
                continue

            outfile.write(f"{chrom}\t{psite_start}\t{psite_end}\t{name}\t{score}\t{strand}\n")
            reads_written += 1

    sys.stderr.write(f"\n=== P-site Correction Statistics for {sample_name} ===\n")
    sys.stderr.write(f"Total reads processed: {reads_processed}\n")
    sys.stderr.write(f"Reads written: {reads_written}\n")
    sys.stderr.write(f"Reads skipped (no offset for length): {reads_skipped_no_offset}\n")
    sys.stderr.write(f"Reads skipped (invalid position): {reads_skipped_invalid_pos}\n")

    if reads_skipped_no_offset > 0:
        sys.stderr.write(f"\nAvailable offsets for read lengths: {sorted(offsets.keys())}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Apply P-site offset correction to deduplicated BED files'
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input deduplicated BED file (6 columns)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output P-site corrected BED file')
    parser.add_argument('-c', '--offsets', required=True,
                        help='CSV file with P-site offsets')
    parser.add_argument('-e', '--experiment', required=True,
                        help='Experiment ID (e.g., GSM2817683)')
    parser.add_argument('-s', '--sample', required=True,
                        help='Sample name for reporting')

    args = parser.parse_args()

    offsets = load_offsets(args.offsets, args.experiment)
    apply_psite_correction_bed(args.input, args.output, offsets, args.sample)

    sys.stderr.write(f"\nP-site correction completed successfully.\n")
    sys.stderr.write(f"Output written to: {args.output}\n")


if __name__ == '__main__':
    main()
