#!/usr/bin/env python3
# Fixed version of hisat2_read_statistics.py for RiboFlow
# Based on original script by Chanhee Park and Daehwan Kim
# Fixed for Python 3 compatibility and empty line handling

import sys, gzip, bz2
from argparse import ArgumentParser

# Constants
COMPRESSION_NON   = 0
COMPRESSION_GZIP  = 1
COMPRESSION_BZIP2 = 2

SEQUENCE_UNKNOWN  = -1
SEQUENCE_FASTA    = 0
SEQUENCE_FASTQ    = 1

FASTA_EXTENSIONS = ["fa", "fasta", "fna"]
FASTQ_EXTENSIONS = ["fq", "fastq"]
MAX_SKIP_LINES = 10000

def parser_FQ(fp):
    skip_line_count = 0
    # Robustly skip empty initial lines
    while skip_line_count < MAX_SKIP_LINES:
        pos = fp.tell()
        line = fp.readline()
        if not line: break 
        
        # Check if line is not empty before indexing
        if line.strip() and line.startswith('@'):
            fp.seek(pos) # Go back to start of header
            break
        skip_line_count += 1
            
    if skip_line_count == MAX_SKIP_LINES:
        raise ValueError("Invalid file format or too many empty lines")

    while True:
        line = fp.readline()
        if not line: return
        
        # Parse standard FASTQ 4-line block
        if not line.startswith('@'): continue
            
        header = line[1:].split()[0]
        seq = fp.readline().strip()
        fp.readline() # +
        fp.readline() # qual
        
        if not seq: return
        yield header, seq

def parse_type(fname):
    compression_type = COMPRESSION_NON
    sequence_type = SEQUENCE_UNKNOWN
    
    if fname.endswith('.gz'):
        compression_type = COMPRESSION_GZIP
        fname = fname[:-3]
    elif fname.endswith('.bz2'):
        compression_type = COMPRESSION_BZIP2
        fname = fname[:-4]
        
    ext = fname.split('.')[-1].lower()
    if ext in FASTA_EXTENSIONS: sequence_type = SEQUENCE_FASTA
    elif ext in FASTQ_EXTENSIONS: sequence_type = SEQUENCE_FASTQ
    
    return sequence_type, compression_type

def reads_stat(read_file, read_count):
    length_map = {}
    
    try:
        seq_type, comp_type = parse_type(read_file)
        
        if comp_type == COMPRESSION_GZIP:
            fp = gzip.open(read_file, 'rt')
        elif comp_type == COMPRESSION_BZIP2:
            fp = bz2.open(read_file, 'rt')
        else:
            fp = open(read_file, 'r')
            
        if seq_type == SEQUENCE_FASTQ:
            stream = parser_FQ(fp)
        else:
            # Fallback for FASTA or unknown (simplified)
            fp.close()
            return

        cnt = 0
        for _, seq in stream:
            l = len(seq)
            length_map[l] = length_map.get(l, 0) + 1
            cnt += 1
            if read_count > 0 and cnt >= read_count: break
            
        fp.close()
        
        if cnt == 0:
            print(f"0 reads checked from {read_file}")
            return

        mn = min(length_map.keys())
        mx = max(length_map.keys())
        
        total_len = sum(k*v for k,v in length_map.items())
        avg = total_len // cnt
        
        # Sort by count (desc), then length
        sorted_map = sorted(length_map.items(), key=lambda x: (x[1], x[0]), reverse=True)
        top_lengths = ",".join(str(k) for k,_ in sorted_map[:10]) # Show top 10
        
        print(f"{cnt} reads; min_len: {mn}; max_len: {mx}; avg_len: {avg}; common_lengths: {top_lengths}")

    except Exception as e:
        sys.stderr.write(f"Error calculating read stats: {str(e)}\n")

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('read_file', nargs='?', help='Input read file')
    parser.add_argument('-n', dest='read_count', type=int, default=10000)
    args = parser.parse_args()
    
    if args.read_file:
        reads_stat(args.read_file, args.read_count)
