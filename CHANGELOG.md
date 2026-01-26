# RELEASE NOTES

## Version 2.0.0 (Current)

### Features
- **Genome Alignment**: Replaced Bowtie2 transcriptome alignment with HISAT2 genome alignment.
- **P-site Offset Correction**:
  - Optional P-site correction using an offset file (`psite_offset` in YAML).
  - Generates P-site corrected BAM and BED files.
  - Updates downstream statistics based on P-site correction.
- **BigWig Generation**:
  - Ribo-Seq: Generates strand-specific BigWig files (requires P-site correction).
  - RNA-Seq: Generates strand-specific BigWig files representing read coverage.
- **Library Strandedness**: Automatic or manual detection (`library_strandedness: "auto"|"forward"|"reverse"`).
- **Genome Read Trimming**: Optional trimming of reads after filtering (`genome_read_trim`).
- **RNA-Seq Processing**: Full support for RNA-Seq data with genome alignment and deduplication options.
- **Statistics**: Comprehensive statistics compilation for genome alignment, including P-site counts.
- **New Scripts**:
  - `apply_psite_offsets.py` / `apply_psite_offsets_bed.py`: Apply P-site offsets to BAM/BED.
  - `detect_strand.py`: Automatically detects library strandedness.
  - `extract_reads_from_dedup_bed.py`: Creates BAM files from deduplicated BEDs.
  - Statistics update scripts: `update_merged_stats_with_counts.py`, `update_merged_stats_with_psite_only.py`, `update_rnaseq_merged_stats.py`, `hisat2_read_statistics_fixed.py`.

### Changes
- **Pipeline Logic**:
  - Default `do_ribo_creation` set to `false`.
  - Output directory structure updated for genome alignment intermediates and results.
  - Ribo file generation is currently disabled/unsupported in this version.
- **Configuration**:
  - New YAML parameters for `psite_offset`, `rnaseq` (with `genome_read_trim`, `library_strandedness`), and updated `alignment_arguments`.
  - Removed RiboR/RiboPy/ribo file support mentions (temporarily).
- **Environment**:
  - Renamed Conda environment to `ribo_genome`.
  - Docker image updated to `danielnguyener/riboflow`.

---

## Version 0.0.1

### Features
- Added UMI support
- Basic UMI extraction and deduplication

---

## Version 0.0.0

### Initial Release
- Initial RiboFlow pipeline release
- Transcriptome alignment support
- Basic deduplication methods
- RIBO file format support
