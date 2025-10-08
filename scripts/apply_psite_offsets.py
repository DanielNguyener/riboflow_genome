#!/usr/bin/env python3

import argparse
import csv
import sys
import pysam


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


def apply_psite_correction(input_bam, output_bam, offsets):
    inbam = pysam.AlignmentFile(input_bam, "rb")
    outbam = pysam.AlignmentFile(output_bam, "wb", template=inbam)

    for read in inbam:
        if read.is_unmapped or read.query_length is None:
            continue

        read_length = read.query_length

        if read_length not in offsets:
            continue

        offset = offsets[read_length]

        if read.is_reverse:
            new_end = read.reference_end - offset
            new_start = new_end - 1
        else:
            new_start = read.reference_start + offset
            new_end = new_start + 1

        if new_start < 0 or new_end <= new_start:
            continue

        new_read = pysam.AlignedSegment()
        new_read.query_name = read.query_name
        new_read.query_sequence = read.query_sequence[0] if read.query_sequence else "N"
        new_read.flag = read.flag
        new_read.reference_id = read.reference_id
        new_read.reference_start = new_start
        new_read.mapping_quality = read.mapping_quality
        new_read.cigar = [(0, 1)]
        new_read.next_reference_id = read.next_reference_id
        new_read.next_reference_start = read.next_reference_start
        new_read.template_length = read.template_length

        if read.query_qualities is not None and len(read.query_qualities) > 0:
            new_read.query_qualities = pysam.qualitystring_to_array(chr(read.query_qualities[0] + 33))

        if read.has_tag('NH'):
            new_read.set_tag('NH', read.get_tag('NH'))
        if read.has_tag('AS'):
            new_read.set_tag('AS', read.get_tag('AS'))
        if read.has_tag('nM'):
            new_read.set_tag('nM', read.get_tag('nM'))

        outbam.write(new_read)

    inbam.close()
    outbam.close()


def main():
    parser = argparse.ArgumentParser(description='Apply P-site offset correction to deduplicated BAM files')

    parser.add_argument('-i', '--input', required=True, help='Input deduplicated BAM file')
    parser.add_argument('-o', '--output', required=True, help='Output P-site corrected BAM file')
    parser.add_argument('-c', '--offsets', required=True, help='CSV file with P-site offsets')
    parser.add_argument('-e', '--experiment', required=True, help='Experiment ID (e.g., GSM2667787)')
    parser.add_argument('-s', '--sample', required=True, help='Sample name for reporting')
    parser.add_argument('--index', action='store_true', help='Create BAM index after correction')

    args = parser.parse_args()

    offsets = load_offsets(args.offsets, args.experiment)
    apply_psite_correction(args.input, args.output, offsets)

    if args.index:
        sorted_bam = args.output.replace('.bam', '.sorted.bam')
        pysam.sort("-o", sorted_bam, args.output)

        import os
        os.replace(sorted_bam, args.output)

        pysam.index(args.output)


if __name__ == '__main__':
    main()
