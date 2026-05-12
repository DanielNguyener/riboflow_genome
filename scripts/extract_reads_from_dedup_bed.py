#!/usr/bin/env python3

import pysam
import sys
import argparse
import os


def extract_one_read_per_bed_coord(bam_file, bed_file, output_bam):
    """
    Extract ONE read per BED coordinate from BAM.
    
    Args:
        bam_file: Input BAM file (quality-filtered, before deduplication)
        bed_file: Deduplicated BED file with unique coordinates
        output_bam: Output BAM file with one read per BED coordinate
    
    Returns:
        Number of reads extracted
    """
    inbam = pysam.AlignmentFile(bam_file, "rb")
    outbam = pysam.AlignmentFile(output_bam, "wb", template=inbam)
    
    bed_coords = {}
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) >= 6:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                strand = fields[5]
                coord_key = (chrom, start, end, strand)
                bed_coords[coord_key] = True
    
    if not bed_coords:
        print("Warning: No valid coordinates found in BED file", file=sys.stderr)
        inbam.close()
        outbam.close()
        return 0
    
    reads_written = 0
    coords_found = 0
    
    for read in inbam:
        if read.is_unmapped:
            continue
        
        chrom = inbam.get_reference_name(read.reference_id)
        start = read.reference_start
        end = read.reference_end
        strand = '-' if read.is_reverse else '+'
        coord_key = (chrom, start, end, strand)
        
        if coord_key in bed_coords and bed_coords[coord_key]:
            outbam.write(read)
            bed_coords[coord_key] = False
            reads_written += 1
            coords_found += 1
    
    coords_not_found = sum(1 for v in bed_coords.values() if v)
    if coords_not_found > 0:
        print(f"Warning: {coords_not_found} coordinates from BED not found in BAM", file=sys.stderr)
    
    inbam.close()
    outbam.close()
    
    print(f"Extracted {reads_written} reads matching {len(bed_coords)} unique BED coordinates", file=sys.stderr)
    
    return reads_written


def main():
    parser = argparse.ArgumentParser(
        description='Extract one read per BED coordinate from BAM file, preserving deduplication'
    )
    parser.add_argument(
        '--bam',
        required=True,
        help='Input BAM file (quality-filtered, before deduplication)'
    )
    parser.add_argument(
        '--bed',
        required=True,
        help='Deduplicated BED file with unique coordinates'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output BAM file with one read per BED coordinate'
    )
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.bam):
        print(f"Error: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.bed):
        print(f"Error: BED file not found: {args.bed}", file=sys.stderr)
        sys.exit(1)
    
    try:
        reads_written = extract_one_read_per_bed_coord(args.bam, args.bed, args.output)
        if reads_written == 0:
            print("Error: No reads extracted", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

