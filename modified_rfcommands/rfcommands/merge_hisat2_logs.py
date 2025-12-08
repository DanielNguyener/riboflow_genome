#!/usr/bin/env python3
"""
Merge multiple hisat2 log files by summing their statistics
"""

import argparse
import re
import sys

def parse_hisat2_log(log_file):
    """Parse hisat2 log file and extract alignment statistics"""
    stats = {
        'total_reads': 0,
        'unaligned': 0,
        'aligned_once': 0,
        'aligned_many': 0
    }

    with open(log_file, 'r') as f:
        content = f.read()

    # Extract total reads
    total_reads_match = re.search(r'(\d+) reads; of these:', content)
    if total_reads_match:
        stats['total_reads'] = int(total_reads_match.group(1))

    # Extract alignment statistics
    unaligned_match = re.search(r'(\d+) \(([\d.]+)%\) aligned 0 times', content)
    if unaligned_match:
        stats['unaligned'] = int(unaligned_match.group(1))

    aligned_once_match = re.search(r'(\d+) \(([\d.]+)%\) aligned exactly 1 time', content)
    if aligned_once_match:
        stats['aligned_once'] = int(aligned_once_match.group(1))

    aligned_many_match = re.search(r'(\d+) \(([\d.]+)%\) aligned >1 times', content)
    if aligned_many_match:
        stats['aligned_many'] = int(aligned_many_match.group(1))

    return stats

def merge_stats(log_files):
    """Merge statistics from multiple log files"""
    merged = {
        'total_reads': 0,
        'unaligned': 0,
        'aligned_once': 0,
        'aligned_many': 0
    }

    for log_file in log_files:
        stats = parse_hisat2_log(log_file)
        for key in merged:
            merged[key] += stats[key]

    return merged

def write_merged_log(stats, output_file):
    """Write merged statistics in hisat2 log format"""
    total_reads = stats['total_reads']
    unaligned = stats['unaligned']
    aligned_once = stats['aligned_once']
    aligned_many = stats['aligned_many']

    # Calculate percentages
    unaligned_pct = (unaligned / total_reads * 100) if total_reads > 0 else 0
    aligned_once_pct = (aligned_once / total_reads * 100) if total_reads > 0 else 0
    aligned_many_pct = (aligned_many / total_reads * 100) if total_reads > 0 else 0
    overall_pct = ((aligned_once + aligned_many) / total_reads * 100) if total_reads > 0 else 0

    with open(output_file, 'w') as f:
        f.write(f"{total_reads} reads; of these:\n")
        f.write(f"  {total_reads} (100.00%) were unpaired; of these:\n")
        f.write(f"    {unaligned} ({unaligned_pct:.2f}%) aligned 0 times\n")
        f.write(f"    {aligned_once} ({aligned_once_pct:.2f}%) aligned exactly 1 time\n")
        f.write(f"    {aligned_many} ({aligned_many_pct:.2f}%) aligned >1 times\n")
        f.write(f"{overall_pct:.2f}% overall alignment rate\n")

def main():
    parser = argparse.ArgumentParser(description='Merge hisat2 log files')
    parser.add_argument('-o', '--output', required=True, help='Output merged log file')
    parser.add_argument('logs', nargs='+', help='Input hisat2 log files to merge')

    args = parser.parse_args()

    try:
        merged_stats = merge_stats(args.logs)
        write_merged_log(merged_stats, args.output)
        print(f"Merged {len(args.logs)} log files into {args.output}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()

