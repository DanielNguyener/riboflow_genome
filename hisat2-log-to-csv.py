#!/usr/bin/env python3
"""
Custom hisat2 log to CSV converter
"""

import argparse
import re
import sys

def parse_hisat2_log(log_file):
    """Parse hisat2 log file and extract alignment statistics"""
    stats = {}
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    total_reads_match = re.search(r'(\d+) reads; of these:', content)
    if total_reads_match:
        stats['total_reads'] = int(total_reads_match.group(1))
    
    unaligned_match = re.search(r'(\d+) \(([\d.]+)%\) aligned 0 times', content)
    if unaligned_match:
        stats['unaligned'] = int(unaligned_match.group(1))
        stats['unaligned_pct'] = float(unaligned_match.group(2))
    
    aligned_once_match = re.search(r'(\d+) \(([\d.]+)%\) aligned exactly 1 time', content)
    if aligned_once_match:
        stats['aligned_once'] = int(aligned_once_match.group(1))
        stats['aligned_once_pct'] = float(aligned_once_match.group(2))
    
    aligned_many_match = re.search(r'(\d+) \(([\d.]+)%\) aligned >1 times', content)
    if aligned_many_match:
        stats['aligned_many'] = int(aligned_many_match.group(1))
        stats['aligned_many_pct'] = float(aligned_many_match.group(2))
    
    if 'aligned_once' in stats and 'aligned_many' in stats:
        stats['total_aligned'] = stats['aligned_once'] + stats['aligned_many']
        
    overall_match = re.search(r'([\d.]+)% overall alignment rate', content)
    if overall_match:
        stats['total_aligned_pct'] = float(overall_match.group(1))
    
    return stats

def write_csv(stats, output_file, sample_name, prefix):
    """Write statistics to CSV file in format compatible with bt2-log-to-csv"""
    
    with open(output_file, 'w') as f:
        f.write(f",{sample_name}\n")
        
        f.write(f"{prefix}_aligned_once,{stats.get('aligned_once', 0)}\n")
        f.write(f"{prefix}_aligned_once_%,{stats.get('aligned_once_pct', 0.0)}\n")
        f.write(f"{prefix}_aligned_many,{stats.get('aligned_many', 0)}\n")
        f.write(f"{prefix}_aligned_many_%,{stats.get('aligned_many_pct', 0.0)}\n")
        f.write(f"{prefix}_total_aligned,{stats.get('total_aligned', 0)}\n")
        f.write(f"{prefix}_total_aligned_%,{stats.get('total_aligned_pct', 0.0)}\n")
        f.write(f"{prefix}_unaligned,{stats.get('unaligned', 0)}\n")
        f.write(f"{prefix}_unaligned_%,{stats.get('unaligned_pct', 0.0)}\n")

def main():
    parser = argparse.ArgumentParser(description='Convert hisat2 alignment log to CSV format')
    parser.add_argument('-l', '--log', required=True, help='hisat2 log file')
    parser.add_argument('-n', '--name', required=True, help='Name of the experiment (column label)')
    parser.add_argument('-p', '--prefix', required=True, help='Prefix for row names')
    parser.add_argument('-o', '--out', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    try:
        stats = parse_hisat2_log(args.log)
        
        write_csv(stats, args.out, args.name, args.prefix)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
