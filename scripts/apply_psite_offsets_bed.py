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


def apply_psite_correction_bed(input_bed, output_bed, offsets):
    try:
        with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
            for line in infile:
                if line.startswith('#') or line.strip() == '':
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue

                try:
                    # BED6 format: chr, start, end, name, score, strand
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    name = fields[3]
                    score = fields[4]
                    strand = fields[5]

                    # Calculate read length from BED coordinates
                    read_length = end - start

                    if read_length in offsets:
                        offset = offsets[read_length]

                        # Apply strand-specific P-site offset (like BAM version)
                        if strand == '-':
                            # For reverse strand: P-site is at the end of the read minus offset
                            psite_pos = end - offset - 1  # -1 because BED is 0-based
                        else:
                            # For forward strand: P-site is at start + offset
                            psite_pos = start + offset

                        # Skip if P-site position is negative
                        if psite_pos < 0:
                            continue

                        # Create P-site BED record (single base)
                        psite_end = psite_pos + 1
                        psite_name = f"{name}_psite"
                        psite_score = score

                        outfile.write(f"{chrom}\t{psite_pos}\t{psite_end}\t{psite_name}\t{psite_score}\t{strand}\n")

                except (ValueError, IndexError):
                    # Skip malformed lines
                    continue

    except FileNotFoundError:
        sys.stderr.write(f"ERROR: Input BED file not found: {input_bed}\n")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description='Apply P-site offset correction to deduplicated BED files')

    parser.add_argument('-i', '--input', required=True, help='Input deduplicated BED file')
    parser.add_argument('-o', '--output', required=True, help='Output P-site corrected BED file')
    parser.add_argument('-c', '--offsets', required=True, help='CSV file with P-site offsets')
    parser.add_argument('-e', '--experiment', required=True, help='Experiment ID (e.g., GSM2667787)')
    parser.add_argument('-s', '--sample', required=True, help='Sample ID (e.g., GSM2667787)')

    args = parser.parse_args()

    # Load offsets for the specified experiment
    offsets = load_offsets(args.offsets, args.experiment)

    # Apply P-site correction to the BED file
    apply_psite_correction_bed(args.input, args.output, offsets)

    sys.stderr.write(f"P-site correction applied successfully for sample {args.sample}\n")
    sys.stderr.write(f"Offsets applied for read lengths: {list(offsets.keys())}\n")


if __name__ == "__main__":
    main()