# -*- coding: utf-8 -*-

from collections import OrderedDict
import pandas as pd

#############################################################################

TOTAL_INDEX = 0
UNAL_INDEX = 2
ONCE_INDEX = 3
MANY_INDEX = 4

def get_reads_from_cutadapt_log(log_file):
    total_reads = -1
    reads_written = -1
    with open(log_file) as input_stream:
      for this_line in input_stream:
          if this_line.startswith("Total reads"):
              total_reads = int( "".join( this_line.split()[-1].split(",") ) )

          if this_line.startswith("Reads written"):
            reads_written = int( "".join( this_line.split()[-2].split(",") ) )
    if total_reads <= 0:
        raise IOError("Couldn't get total reads from the cutadapt log {}".\
                        format(log_file))
    if reads_written <= 0:
        raise IOError("Couldn't get reads_written from the cutadapt log {}".\
                        format(log_file))
    return (total_reads, reads_written)


def read_bowtie2_log(log_file):
    log_lines = []
    with open(log_file) as input_stream:
        for this_line in input_stream:
            # Ignore all warnings etc
            # The actual log starts with a numeric entry
            # The warnings and other messages start with a non-numeric character
            if len(this_line) < 1:
                continue
            if this_line[0].isalpha():
                continue
            log_lines.append(this_line)
    if len(log_lines) != 6:
        raise IOError("The file {} has to contain exactly 6 lines".format(log_file) )
    return list( map( lambda this_line: int( this_line.split()[0] ) ,log_lines[:5] ) )

def read_hisat2_log(log_file):
    """Parse hisat2 alignment log file and extract statistics in bowtie2 format"""
    total_reads = 0
    unaligned = 0
    once_aligned = 0
    multi_aligned = 0

    with open(log_file) as input_stream:
        for this_line in input_stream:
            line = this_line.strip()
            # Parse hisat2 output format
            if "reads; of these:" in line:
                total_reads = int(line.split()[0])
            elif "aligned 0 times" in line:
                unaligned = int(line.split()[0])
            elif "aligned exactly 1 time" in line:
                once_aligned = int(line.split()[0])
            elif "aligned >1 times" in line:
                multi_aligned = int(line.split()[0])

    # Return in same format as bowtie2: [total, ?, unaligned, once, multi]
    return [total_reads, 0, unaligned, once_aligned, multi_aligned]

def detect_log_type(log_file):
    """Detect if log file is from bowtie2 or hisat2"""
    with open(log_file) as input_stream:
        content = input_stream.read()
        if "overall alignment rate" in content and "aligned exactly 1 time" in content:
            return "hisat2"
        else:
            return "bowtie2"

def get_count_from_qpass_file(qpass_file):
    with open(qpass_file) as input_stream:
        qpassing_align_count = input_stream.readlines()[0].strip()
    return int(qpassing_align_count)

def get_count_from_dedup(dedup_file):
    with open(dedup_file, 'r') as input_stream:
        this_line = input_stream.readline().strip().split()[0]
    return int(this_line)


def get_overall_statistics(cutadapt_log, filter_log, transcriptome_log,
                           qpass_count_file, dedup_file, label_prefix="transcriptome", psite_file=None):
    overall_statistics = OrderedDict()
    overall_statistics["total_reads"], overall_statistics["clipped_reads"] = \
          get_reads_from_cutadapt_log(cutadapt_log)

    filter_stats = read_bowtie2_log(filter_log)
    overall_statistics["filtered_out"] = filter_stats[ONCE_INDEX] + filter_stats[MANY_INDEX]
    overall_statistics["filter_kept"]  = filter_stats[UNAL_INDEX]

    # Detect log type and use appropriate reader
    log_type = detect_log_type(transcriptome_log)
    if log_type == "hisat2":
        transcriptome_stats = read_hisat2_log(transcriptome_log)
    else:
        transcriptome_stats = read_bowtie2_log(transcriptome_log)
    overall_statistics[f"{label_prefix}_aligned_once"] = transcriptome_stats[ONCE_INDEX]

    overall_statistics[f"{label_prefix}_aligned_many"] = transcriptome_stats[MANY_INDEX]

    overall_statistics[f"{label_prefix}_total_aligned"] = \
         overall_statistics[f"{label_prefix}_aligned_once"] + \
         overall_statistics[f"{label_prefix}_aligned_many"]

    overall_statistics[f"{label_prefix}_unaligned"] = transcriptome_stats[UNAL_INDEX]

    overall_statistics[f"{label_prefix}_qpass_aligned_reads"] = get_count_from_qpass_file(qpass_count_file)

    overall_statistics[f"{label_prefix}_after_dedup"] = get_count_from_dedup(dedup_file)

    if psite_file:
        overall_statistics[f"{label_prefix}_psite_hits"] = get_count_from_dedup(psite_file)

    """
    genome_stats = read_bowtie2_log(genome_log)
    overall_statistics["genome_aligned_once"] = genome_stats[ONCE_INDEX]

    overall_statistics["genome_aligned_many"] = genome_stats[MANY_INDEX]

    overall_statistics["genome_total_aligned"] = \
         overall_statistics["genome_aligned_once"] + \
         overall_statistics["genome_aligned_many"]

    overall_statistics["genome_unaligned"] = genome_stats[UNAL_INDEX]
    """
    return overall_statistics

def print_overall_statistics(stats_dict, sample_name, csv_file):
	stats_df = pd.DataFrame.from_dict(stats_dict,
                                      orient  = 'index',
                                      columns = [sample_name])
	stats_df.to_csv(csv_file)


def compile(out,     cutadapt, filter, trans,
            quality,   dedup,  name, label_prefix="transcriptome", psite=None):

    overall_statistics = get_overall_statistics(
                                cutadapt_log      = cutadapt,
                                filter_log        = filter,
                                transcriptome_log = trans,
                                qpass_count_file  = quality,
                                dedup_file        = dedup,
                                label_prefix      = label_prefix,
                                psite_file        = psite)

    print_overall_statistics(stats_dict  = overall_statistics,
    	                     sample_name = name,
    	                     csv_file    = out)