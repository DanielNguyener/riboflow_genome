# RELEASE NOTES

## Version 2.0.0 (Current)

### Features
- Genome alignment using HISAT2
- P-site offset correction for BED and BAM files
- Strand-specific bigWig generation
- RNA-seq processing with genome alignment
- Statistics compilation with genome alignment support
- New commands: `rfc hisat2-log-to-csv`, `rfc merge-hisat2-logs`
- New scripts: `apply_psite_offsets_bed.py`, `apply_psite_offsets.py`, `extract_reads_from_dedup_bed.py`, `update_merged_stats_with_counts.py`, `update_merged_stats_with_psite_only.py`, `update_rnaseq_merged_stats.py`
- New YAML parameters for P-site offsets and RNA-seq configuration
- New conda environment: `ribo_bigwig`

### Changes
- Replaced transcriptome alignment with genome alignment
- Requires HISAT2 genome index instead of Bowtie2 transcriptome index
- Statistics format changed (genome vs transcriptome labels)
- Output directory structure changed
- `do_ribo_creation` default changed to `false`

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
