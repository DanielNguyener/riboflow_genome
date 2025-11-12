# RiboFlow_Genome: Changes from Reference Version

This document lists the differences between RiboFlow_Genome and the reference RiboFlow version (`rf_test_run/riboflow/RiboFlow.groovy`).

## Version Comparison

| Metric | Reference Version | Genome Version |
|--------|------------------|----------------|
| Total Lines | 2,256 | 3,253 |
| Processes | ~25 | 42 |
| Alignment Method | Transcriptome (Bowtie2) | Genome (HISAT2) |
| P-site Correction | No | Yes |
| BigWig Generation | No | Yes |
| RNA-seq Support | Limited | Full |

## Changes

### 1. Alignment Method: Transcriptome → Genome

**Reference Version:**
- Uses Bowtie2 for transcriptome alignment
- Aligns reads to transcriptome reference

**Genome Version:**
- Uses HISAT2 for genome alignment
- Aligns reads to genome reference with splice site awareness
- Transcriptome alignment removed
- Requires HISAT2 genome index instead of Bowtie2 transcriptome index

### 2. P-site Offset Correction

**Added:**
- `apply_psite_offsets_bed.py`: P-site correction for BED files
- `apply_psite_offsets.py`: P-site correction for BAM files
- P-site offset CSV file support (`outputHuman_OFFSETS.csv`)
- Sample-to-experiment mapping in YAML configuration

**Processes Added:**
- `apply_psite_correction_umi_individual`: P-site correction for UMI deduplication
- `apply_psite_correction_none_individual`: P-site correction for 'none' deduplication
- `merge_psite_corrected_bam`: Merge P-site corrected BAM files
- `merge_psite_corrected_bed`: Merge P-site corrected BED files
- `count_individual_psite_bed_none`: Count P-site reads

**Output:**
- P-site corrected BAM files: 1bp reads at P-site positions
- P-site corrected BED files: 1bp intervals at P-site positions

### 3. Strand-Specific BigWig Generation

**Added:**
- `genome_create_strand_specific_bigwigs`: Ribo-seq bigWig generation
- `rnaseq_create_dedup_bigwigs`: RNA-seq bigWig (with deduplication)
- `rnaseq_create_nodedup_bigwigs`: RNA-seq bigWig (without deduplication)
- Separate conda environment: `ribo_bigwig` (for deeptools)

**Output Files:**
- Ribo-seq: `{sample}.psite.plus.bigWig`, `{sample}.psite.minus.bigWig`
- RNA-seq (dedup): `{experiment}.rnaseq.dedup.plus.bigWig`, `{experiment}.rnaseq.dedup.minus.bigWig`
- RNA-seq (no dedup): `{sample}.rnaseq.nodedup.plus.bigWig`, `{sample}.rnaseq.nodedup.minus.bigWig`

**Library Strandedness:**
- Configurable via `library_strandedness` parameter ("reverse" or "forward")

### 4. RNA-seq Support

**Reference Version:**
- Limited RNA-seq processing
- Basic filtering and alignment

**Genome Version:**
- RNA-seq pipeline with genome alignment
- Position-based deduplication support
- Strand-specific bigWig generation
- Separate statistics tracking
- Optional read trimming before alignment

**New RNA-seq Processes:**
- `rnaseq_raw_fastqc`: Quality control
- `rnaseq_clip`: Adapter clipping
- `rnaseq_filter`: rRNA filtering
- `rnaseq_genome_read_trim`: Optional read trimming
- `rnaseq_genome_alignment`: Genome alignment
- `rnaseq_genome_quality_filter`: Quality filtering
- `rnaseq_genome_deduplicate_position`: Position-based deduplication
- `rnaseq_create_dedup_bigwigs`: BigWig generation (with dedup)
- `rnaseq_create_nodedup_bigwigs`: BigWig generation (no dedup)
- RNA-seq statistics compilation processes

### 5. Statistics System

**Reference Version:**
- Statistics for transcriptome alignment
- Basic step statistics

**Genome Version:**
- Statistics for genome alignment
- `--label-prefix` option added
- Separate individual and merged statistics
- P-site count tracking
- RNA-seq statistics tracking

**New Statistics Commands:**
- `rfc hisat2-log-to-csv`: Convert HISAT2 logs to CSV
- `rfc merge-hisat2-logs`: Merge HISAT2 logs
- `rfc compile-step-stats` with `--label-prefix` option
- `rfc stats-percentage` with `--label-prefix` option

**Statistics Update Scripts:**
- `update_merged_stats_with_counts.py`: Update with dedup and P-site counts
- `update_merged_stats_with_psite_only.py`: Update with P-site count only
- `update_rnaseq_merged_stats.py`: Update RNA-seq statistics

### 6. Modified RiboFlow Commands (rfc)

**New Commands:**
- `rfc hisat2-log-to-csv`: HISAT2 log conversion
- `rfc merge-hisat2-logs`: HISAT2 log merging

**Modified Commands:**
- `rfc compile-step-stats`: Added `--label-prefix` option
- `rfc stats-percentage`: Added `--label-prefix` option

**Unused Commands:**
- `rfc bt2-log-to-csv`: Available but not used

### 7. Configuration Changes

**New YAML Parameters:**
```yaml
# P-site offset configuration
psite_offset:
  offset_file: "./outputHuman_OFFSETS.csv"
  sample_matching:
    experiment_1:
      lanes:
        1: "GSM4114767"
        2: "GSM4114768"

# RNA-seq configuration
rnaseq:
  dedup_method: "position"  # or "none"
  library_strandedness: "reverse"  # or "forward"
  genome_read_trim:
    enabled: true
    length: 35

# Disable ribo creation for genome mode
do_ribo_creation: false
```

**Removed/Deprecated:**
- Transcriptome alignment parameters (deprecated)
- `do_ribo_creation` default changed to `false`

### 8. Directory Structure Changes

**New Directories:**
- `intermediates_*/genome_alignment/`: Genome alignment outputs
- `intermediates_*/deduplication/bigwigs/`: BigWig files
- `intermediates_*/rnaseq/`: RNA-seq specific intermediates

**Output Structure:**
```
intermediates_umi/
├── genome_alignment/
│   ├── individual/
│   └── merged/
├── deduplication/
│   ├── merged/
│   └── bigwigs/merged/
└── rnaseq/
    ├── genome_alignment/
    ├── deduplication/
    └── bigwigs/merged/
```

### 9. Function Additions

**New Functions:**
- `get_rnaseq_storedir(output_type)`: RNA-seq storage directory
- `get_rnaseq_publishdir(output_type)`: RNA-seq publish directory
- `get_gsm_id_for_lane(sample, index, file_exists_check)`: Get experiment ID for P-site lookup
- `psite_offset_file_exists`: Check for P-site offset file

**Modified Functions:**
- `get_storedir(output_type, is_rnaseq)`: Added RNA-seq support
- `get_publishdir(output_type, is_rnaseq)`: Added RNA-seq support

**Unused Functions:**
- `build_samtools_sort_cmd`: Defined but not used
- `build_samtools_index_cmd`: Defined but not used
- `build_cutadapt_cmd`: Defined but not used
- `build_fastqc_cmd`: Defined but not used

### 10. Channel Architecture

**New Channels:**
- `GENOME_*`: Genome alignment channels
- `GENOME_PSITE_*`: P-site correction channels
- `GENOME_BIGWIG_*`: BigWig generation channels
- `RNASEQ_*`: RNA-seq specific channels
- `RNASEQ_GENOME_*`: RNA-seq genome alignment channels

**Channel Patterns:**
- Channel splitting for parallel processing
- Conditional channel routing based on deduplication method
- Separate channels for individual vs merged processing

## Backward Compatibility

**Not Backward Compatible:**
- Transcriptome alignment removed
- Statistics format changed (genome vs transcriptome labels)
- Output directory structure changed

**Maintained Compatibility:**
- YAML configuration structure (with additions)
- Deduplication methods (position, umi_tools, none)
- Basic pipeline flow (clip → filter → align → dedup)

## References

- Reference Version: `/root/Cenik/Project_2/vg/processing/2_riboflow/rf_test_run/riboflow/RiboFlow.groovy`
- Genome Version: `/root/Cenik/Project_2/vg/processing/2_riboflow/riboflow/RiboFlow.groovy`
- Code Audit: See `docs/CODE_AUDIT.md`
- Scripts Documentation: See `scripts/README.md`
- Commands Documentation: See `modified_rfcommands/README.md`
