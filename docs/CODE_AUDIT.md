# Code Audit Report

This document lists all processes, channels, functions, commands, and scripts in RiboFlow_Genome and their usage status.

## Summary

- **Total Processes**: 42 processes defined
- **Total Functions**: 11 functions defined
- **Commands in `/modified_rfcommands`**: 8 commands available, 7 used, 1 unused
- **Scripts in `/scripts`**: 6 scripts, all used

## Commands Audit (`/modified_rfcommands`)

### Used Commands

1. **`rfc hisat2-log-to-csv`**
   - Converts HISAT2 alignment logs to CSV format
   - Used at lines: 435, 683, 2550, 2734

2. **`rfc merge-hisat2-logs`**
   - Merges multiple HISAT2 log files
   - Used at line: 682

3. **`rfc dedup`**
   - Position-based deduplication of BED files
   - Used at lines: 770, 2807

4. **`rfc compile-step-stats`**
   - Compiles statistics from multiple pipeline steps
   - Used at lines: 1953, 2888

5. **`rfc stats-percentage`**
   - Calculates percentage statistics from raw counts
   - Used at lines: 2006, 2168, 3109, 3198

6. **`rfc sum-stats`**
   - Sums statistics across multiple lanes/samples
   - Used at lines: 2099, 2111, 2123, 2133, 3153, 3163

7. **`rfc merge overall-stats`**
   - Merges overall statistics tables
   - Used at lines: 2003, 2165, 3106, 3195

### Unused Commands

1. **`rfc bt2-log-to-csv`**
   - Status: Available but not used
   - Purpose: Convert Bowtie2 alignment logs to CSV format

## Scripts Audit (`/scripts`)

All 6 scripts are used:

1. **`apply_psite_offsets_bed.py`**
   - Applies P-site offsets to BED files
   - Used at lines: 1108, 1210

2. **`apply_psite_offsets.py`**
   - Applies P-site offsets to BAM files
   - Used at lines: 1285, 1331

3. **`extract_reads_from_dedup_bed.py`**
   - Extracts one read per BED coordinate from BAM files
   - Used at lines: 1484, 2986

4. **`update_merged_stats_with_counts.py`**
   - Updates merged statistics with dedup and P-site counts
   - Used at lines: 2102, 2114

5. **`update_merged_stats_with_psite_only.py`**
   - Updates merged statistics with P-site count only
   - Used at line: 2126

6. **`update_rnaseq_merged_stats.py`**
   - Updates RNA-seq merged statistics with dedup count
   - Used at line: 3156

## Functions Audit

1. **`get_storedir(output_type, is_rnaseq = false)`**
   - Returns storage directory path for intermediates
   - Used throughout pipeline

2. **`get_rnaseq_storedir(output_type)`**
   - Returns storage directory path for RNA-seq intermediates
   - Used in RNA-seq processes

3. **`get_rnaseq_publishdir(output_type)`**
   - Returns publish directory path for RNA-seq outputs
   - Used in RNA-seq processes

4. **`get_publishdir(output_type, is_rnaseq = false)`**
   - Returns publish directory path for outputs
   - Used throughout pipeline

5. **`get_dedup_method(String dedup_arg, String dedup_old)`**
   - Parses and validates deduplication method
   - Used at line 85-86

6. **`build_samtools_sort_cmd(...)`**
   - Status: Defined but not used

7. **`build_samtools_index_cmd(...)`**
   - Status: Defined but not used

8. **`build_cutadapt_cmd(...)`**
   - Status: Defined but not used

9. **`build_fastqc_cmd(...)`**
   - Status: Defined but not used

10. **`get_gsm_id_for_lane(sample, index, file_exists_check = null)`**
    - Returns GSM/experiment ID for P-site offset lookup
    - Used in P-site correction processes

11. **`psite_offset_file_exists`**
    - Checks if P-site offset file exists
    - Used in conditional P-site correction processes

## Processes Audit

All 42 processes are used (conditional execution based on parameters):

### Ribo-seq Processes
- `write_fastq_correspondence`
- `raw_fastqc`
- `clip`
- `extract_umi_via_umi_tools` (conditional)
- `clipped_fastqc`
- `filter`
- `genome_alignment`
- `genome_quality_filter`
- `merge_genome_bam_post_qpass`
- `individual_genome_bam_to_bed`
- `add_sample_index_col_to_genome_bed`
- `merge_genome_alignment`
- `genome_bam_to_bed`
- `genome_deduplicate_position`
- `genome_deduplicate_umi_tools`
- `genome_umi_dedup_bam_to_bed`
- `apply_psite_correction_umi_individual`
- `apply_psite_correction_none_individual`
- `merge_psite_corrected_bam`
- `count_individual_psite_bed_none`
- `merge_psite_corrected_bed`
- `genome_create_strand_specific_bigwigs`
- `individual_genome_alignment_stats`
- `combine_individual_genome_alignment_stats`
- `sum_individual_genome_alignment_stats`
- `combine_merged_genome_alignment_stats`
- `publish_stats`

### RNA-seq Processes
- `rnaseq_raw_fastqc`
- `rnaseq_clip`
- `rnaseq_filter`
- `rnaseq_genome_read_trim` (conditional)
- `rnaseq_genome_alignment`
- `rnaseq_genome_quality_filter`
- `rnaseq_merge_genome_bam_post_qpass`
- `rnaseq_individual_genome_bam_to_bed`
- `rnaseq_add_sample_index_col_to_genome_bed`
- `rnaseq_merge_genome_alignment`
- `rnaseq_genome_bam_to_bed`
- `rnaseq_genome_deduplicate_position`
- `rnaseq_create_dedup_bigwigs`
- `rnaseq_create_nodedup_bigwigs`
- `rnaseq_individual_genome_alignment_stats`
- `combine_individual_rnaseq_genome_alignment_stats`
- `sum_individual_rnaseq_genome_alignment_stats`
- `combine_merged_rnaseq_genome_alignment_stats`

### Ribo File Creation (Conditional)
- `create_ribo` (when `do_ribo_creation = true`)
- `merge_ribos` (when `do_ribo_creation = true`)

## Channels Audit

All channels are used as part of the data flow. Channel naming patterns:

- Input channels: `INPUT_SAMPLES_*`
- Processing channels: `*_OUT`, `*_LOG`, `*_BAM`, `*_BED`
- Merge channels: `*_GROUPED`, `*_MERGED`
- Output channels: `*_FINAL`, `*_FOR_PUBLISH`
