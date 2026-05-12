# RiboFlow Scripts Documentation

This directory contains Python scripts used by the RiboFlow_Genome pipeline for specialized processing tasks, including P-site correction, read extraction, and statistics updates.

## Overview

All scripts in this directory are actively used by the RiboFlow.groovy pipeline. They handle:
- P-site offset correction for BED and BAM files
- Read extraction from deduplicated BED coordinates
- Statistics CSV updates with deduplication and P-site counts

## Scripts

### 1. `apply_psite_offsets_bed.py`

Applies P-site offset correction to deduplicated BED files.

**Purpose:**
Converts read positions in BED files to P-site positions based on experiment-specific offsets. Used for UMI and position-based deduplication methods.

**Usage:**
```bash
python3 apply_psite_offsets_bed.py \
    -i <input.bed> \
    -o <output.bed> \
    -c <offsets.csv> \
    -e <experiment_id> \
    -s <sample_id>
```

**Arguments:**
- `-i, --input` (required): Input deduplicated BED file
- `-o, --output` (required): Output P-site corrected BED file
- `-c, --offsets` (required): CSV file containing P-site offsets (e.g., `outputHuman_OFFSETS.csv`)
- `-e, --experiment` (required): Experiment ID from the offset CSV (e.g., `GSM4114767`)
- `-s, --sample` (required): Sample ID for reporting

**Input Format:**
- BED file with 6+ columns: `chrom`, `start`, `end`, `name`, `score`, `strand`
- Offset CSV with columns: `Experiment`, `Read Length`, `P-site Offset`

**Output Format:**
- BED file with 1bp intervals at P-site positions
- Each entry converted from full read to P-site position
- Name field appended with `_psite` suffix

**P-site Calculation:**
- **Forward strand (`+`)**: `psite_pos = start + offset`
- **Reverse strand (`-`)**: `psite_pos = end - offset - 1`

**Used in Pipeline:**
- Line 1252: P-site correction for UMI deduplication (individual lanes)
- Line 1354: P-site correction for position deduplication (individual lanes)

**Dependencies:**
- Python 3.7+
- Standard library: `argparse`, `csv`, `sys`

**Error Handling:**
- Exits with error if offset file not found
- Exits with error if experiment ID not found in offset file
- Skips invalid BED entries (continues processing)

---

### 2. `apply_psite_offsets.py`

Applies P-site offset correction to deduplicated BAM files.

**Purpose:**
Converts BAM reads to 1bp reads positioned at the P-site. Used for the 'none' deduplication method where P-site correction is applied directly to BAM files.

**Usage:**
```bash
python3 apply_psite_offsets.py \
    -i <input.bam> \
    -o <output.bam> \
    -c <offsets.csv> \
    -e <experiment_id> \
    -s <sample_name> \
    [--index]
```

**Arguments:**
- `-i, --input` (required): Input deduplicated BAM file
- `-o, --output` (required): Output P-site corrected BAM file
- `-c, --offsets` (required): CSV file containing P-site offsets
- `-e, --experiment` (required): Experiment ID from the offset CSV
- `-s, --sample` (required): Sample name for reporting
- `--index` (optional): Create BAM index after correction

**Input Format:**
- BAM file with aligned reads
- Offset CSV with columns: `Experiment`, `Read Length`, `P-site Offset`

**Output Format:**
- BAM file with 1bp reads at P-site positions
- Each read converted to single base pair at P-site
- Read sequence set to single nucleotide
- CIGAR string set to `1M` (1 match)

**P-site Calculation:**
- **Forward strand**: `new_start = reference_start + offset`, `new_end = new_start + 1`
- **Reverse strand**: `new_end = reference_end - offset`, `new_start = new_end - 1`

**Used in Pipeline:**
- Line 1429: P-site correction for UMI deduplication method (individual lanes)
- Line 1475: P-site correction for 'none' deduplication method (individual lanes)
- Line 1666: P-site correction for position deduplication method (merged, after read extraction)

**Dependencies:**
- Python 3.7+
- `pysam` (for BAM file handling)
- Standard library: `argparse`, `csv`, `sys`

**Error Handling:**
- Exits with error if offset file not found
- Exits with error if experiment ID not found
- Skips unmapped reads and reads without valid offsets
- Skips reads where P-site position would be invalid

---

### 3. `extract_reads_from_dedup_bed.py`

Extracts one representative read per BED coordinate from a BAM file.

**Purpose:**
After deduplication, extracts one read per unique coordinate from the original quality-filtered BAM file. Used for generating bigWig files from deduplicated data.

**Usage:**
```bash
python3 extract_reads_from_dedup_bed.py \
    -b <input.bam> \
    -d <dedup.bed> \
    -o <output.bam>
```

**Arguments:**
- `-b, --bam` (required): Input BAM file (quality-filtered, before deduplication)
- `-d, --dedup-bed` (required): Deduplicated BED file with unique coordinates
- `-o, --output` (required): Output BAM file with one read per BED coordinate

**Input Format:**
- BAM file: Quality-filtered aligned reads
- BED file: 6-column format with deduplicated coordinates

**Output Format:**
- BAM file containing one read per unique BED coordinate
- Reads match the coordinates in the deduplicated BED file

**Algorithm:**
1. Load all coordinates from deduplicated BED file
2. For each read in BAM file:
   - Check if read matches a BED coordinate (chrom, start, end, strand)
   - If match found and no read extracted yet for this coordinate, write read
   - Mark coordinate as processed

**Used in Pipeline:**
- Line 1659: Extract reads for bigWig generation (Ribo-seq, position deduplication, merged)
- Line 1769: Extract reads for bigWig generation (Ribo-seq, position deduplication, no P-site)
- Line 3132: Extract reads for bigWig generation (RNA-seq, position deduplication)

**Dependencies:**
- Python 3.7+
- `pysam` (for BAM file handling)
- Standard library: `argparse`, `sys`, `os`

**Error Handling:**
- Warns if no valid coordinates found in BED file
- Skips invalid BED entries
- Reports number of reads extracted

---

### 4. `update_merged_stats_with_counts.py`

Updates merged statistics CSV with deduplication and P-site counts.

**Purpose:**
Updates the merged statistics CSV file with actual counts from merged deduplication and P-site correction steps. Replaces estimated counts with actual merged counts.

**Usage:**
```bash
python3 update_merged_stats_with_counts.py \
    --dedup-count-file <dedup.count> \
    [--psite-count-file <psite.count>] \
    --input-csv <input.csv> \
    --output-csv <output.csv>
```

**Arguments:**
- `--dedup-count-file` (required): Path to merged deduplication count file (contains single number)
- `--psite-count-file` (optional): Path to merged P-site count file (contains single number)
- `--input-csv` (required): Input temporary CSV file with estimated counts
- `--output-csv` (required): Output CSV file with updated counts

**Input Format:**
- Count files: Plain text files containing a single integer
- CSV file: Statistics CSV with columns including `genome_after_dedup` and optionally `genome_after_psite`

**Output Format:**
- CSV file with updated `genome_after_dedup` column
- If P-site count provided, updates `genome_after_psite` column
- All other columns preserved

**Used in Pipeline:**
- Line 2256: Update merged stats with dedup and P-site counts (UMI deduplication with P-site)
- Line 2268: Update merged stats with dedup count only (UMI/position deduplication without P-site)

**Dependencies:**
- Python 3.7+
- Standard library: `argparse`, `csv`, `sys`

**Error Handling:**
- Exits with error if dedup count file not found
- Gracefully handles missing P-site count file (if optional)
- Exits with error if input CSV not found
- Preserves all CSV columns and structure

---

### 5. `update_merged_stats_with_psite_only.py`

Updates merged statistics CSV with P-site count only.

**Purpose:**
Updates the merged statistics CSV file with P-site count when deduplication is not performed (dedup_method = 'none'). Only updates the P-site count column.

**Usage:**
```bash
python3 update_merged_stats_with_psite_only.py \
    --psite-count-file <psite.count> \
    --input-csv <input.csv> \
    --output-csv <output.csv>
```

**Arguments:**
- `--psite-count-file` (required): Path to merged P-site count file (contains single number)
- `--input-csv` (required): Input temporary CSV file
- `--output-csv` (required): Output CSV file with updated P-site count

**Input Format:**
- Count file: Plain text file containing a single integer
- CSV file: Statistics CSV with `genome_after_psite` column

**Output Format:**
- CSV file with updated `genome_after_psite` column
- All other columns preserved

**Used in Pipeline:**
- Line 2280: Update merged stats with P-site count only (no deduplication method)

**Dependencies:**
- Python 3.7+
- Standard library: `argparse`, `csv`, `sys`

**Error Handling:**
- Exits with error if P-site count file not found
- Exits with error if input CSV not found
- Preserves all CSV columns and structure

---

### 6. `update_rnaseq_merged_stats.py`

Updates RNA-seq merged statistics CSV with deduplication count.

**Purpose:**
Updates the RNA-seq merged statistics CSV file with actual deduplication count from merged RNA-seq BAM files. Similar to `update_merged_stats_with_counts.py` but specifically for RNA-seq data.

**Usage:**
```bash
python3 update_rnaseq_merged_stats.py \
    --dedup-count-file <dedup.count> \
    --input-csv <input.csv> \
    --output-csv <output.csv>
```

**Arguments:**
- `--dedup-count-file` (required): Path to merged RNA-seq deduplication count file
- `--input-csv` (required): Input temporary CSV file with estimated counts
- `--output-csv` (required): Output CSV file with updated counts

**Input Format:**
- Count file: Plain text file containing a single integer
- CSV file: RNA-seq statistics CSV with `genome_after_dedup` column

**Output Format:**
- CSV file with updated `genome_after_dedup` column
- All other columns preserved

**Used in Pipeline:**
- Line 3303: Update RNA-seq merged stats with dedup count (position deduplication)

**Dependencies:**
- Python 3.7+
- Standard library: `argparse`, `csv`, `sys`

**Error Handling:**
- Exits with error if dedup count file not found
- Exits with error if input CSV not found
- Preserves all CSV columns and structure

---

## Common Patterns

### P-site Offset File Format

All P-site correction scripts expect a CSV file with the following format:

```csv
Experiment,Read Length,P-site Offset
GSM4114767,28,12
GSM4114767,29,12
GSM4114767,30,13
...
```

- **Experiment**: Experiment ID (e.g., GSM accession number)
- **Read Length**: Read length in nucleotides
- **P-site Offset**: Offset from 5' end (forward) or 3' end (reverse) to P-site

### Count File Format

Statistics update scripts expect count files to be plain text files containing a single integer:

```
1234567
```

No headers, no formatting, just the count as a single number.

### CSV Statistics Format

Statistics CSV files follow a standard format with columns for each processing step:

- `total_reads`: Total input reads
- `clipped_reads`: Reads after adapter clipping
- `filtered_out`: Reads removed by rRNA filtering
- `filter_kept`: Reads passing rRNA filtering
- `genome_aligned_once`: Reads aligned exactly once
- `genome_aligned_many`: Reads aligned multiple times
- `genome_total_aligned`: Total aligned reads
- `genome_unaligned`: Unaligned reads
- `genome_qpass_aligned_reads`: Quality-filtered aligned reads
- `genome_after_dedup`: Reads after deduplication
- `genome_after_psite`: Reads after P-site correction

## Dependencies

All scripts require:
- **Python**: 3.7 or higher
- **Standard Library**: `argparse`, `csv`, `sys`, `os`

Additional dependencies:
- **Pysam**: Required for `apply_psite_offsets.py` and `extract_reads_from_dedup_bed.py`
  - Install: `pip install pysam` or `conda install pysam`

## Integration with Pipeline

All scripts are automatically called by the RiboFlow.groovy pipeline. They are invoked using:

```groovy
python3 ${workflow.projectDir}/scripts/<script_name>.py <arguments>
```

Scripts are executed within the `ribo_genome` conda environment, which includes all required dependencies.

## Error Handling

All scripts follow consistent error handling patterns:
- Print error messages to `stderr`
- Exit with non-zero status codes on failure
- Provide clear error messages for common issues (file not found, missing columns, etc.)
- Continue processing when possible (e.g., skip invalid entries)

## Testing

To test scripts individually:

```bash
# Activate conda environment
conda activate ribo_genome

# Test P-site correction
python3 scripts/apply_psite_offsets_bed.py --help

# Test read extraction
python3 scripts/extract_reads_from_dedup_bed.py --help

# Test statistics update
python3 scripts/update_merged_stats_with_counts.py --help
```

## Troubleshooting

### Import Errors

If you encounter `ModuleNotFoundError`:

```bash
# Ensure conda environment is activated
conda activate ribo_genome

# Verify pysam installation
python -c "import pysam; print(pysam.__version__)"

# Reinstall if needed
pip install pysam
```

### File Not Found Errors

- Verify input file paths are correct
- Check that files exist and are readable
- Ensure paths are absolute or relative to working directory

### Invalid Format Errors

- Verify BED files have at least 6 columns
- Check that CSV files have required columns
- Ensure count files contain only integers

## Version Compatibility

- **Python**: 3.7+
- **Pysam**: 0.15.0+ (for BAM handling scripts)
- **Nextflow**: 19.04.1+ (for pipeline integration)

## License

See LICENSE file in the main repository.

