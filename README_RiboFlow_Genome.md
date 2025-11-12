# RiboFlow_Genome

A comprehensive Nextflow pipeline for processing ribosome profiling (Ribo-seq) and RNA-seq data with genome alignment, UMI deduplication, P-site correction, and strand-specific bigWig generation.

## Table of Contents

- [Quick Start](#quick-start)
- [Required Files](#required-files)
- [Installation & Dependencies](#installation--dependencies)
- [YAML Configuration](#yaml-configuration)
- [P-site Offset Correction](#p-site-offset-correction-ribo-seq)
- [Strand-Specific BigWig Generation](#strand-specific-bigwig-file-generation)
- [Directory Structure](#directory-structure-overview)
- [Notes](#notes)
- [Version Information](#version-information)

---

## Quick Start

Get up and running with RiboFlow_Genome in three simple steps:

### Environment Setup (First time only)
```bash
# Install conda environments
conda env create -f environment.yaml
conda env create -f ribo_bigwig_environment.yaml
```

### Activate Environment
```bash
conda activate ribo_genome
```

### Run Pipeline
```bash
# Basic run
nextflow run RiboFlow.groovy -params-file project_umi.yaml

# With custom configuration (optional)
nextflow run RiboFlow.groovy -params-file project_umi.yaml -c configs/local.config
```

Make sure your configuration file is properly configured before running the pipeline. See [YAML Configuration](#yaml-configuration) for details.

---

## Required Files

- Genome index files (included in `./rf_sample_data/genome/`)


---
## Installation & Dependencies

### Conda Environments

#### Main Environment: `ribo_genome`
Created from `environment.yaml` - contains all core pipeline tools:
```bash
conda env create -f environment.yaml
```

This environment includes:
- `python=3.7`, `nextflow=19.04.1`, `openjdk=8`
- Alignment tools: `hisat2`, `bowtie2`, `samtools>=1.4`
- Processing tools: `cutadapt`, `bedtools`, `fastqc`
- Python packages: `numpy`, `pandas`, `matplotlib`, `pysam`
- Pip packages: `umi_tools`, `ribopy`, `scipy`, and custom `rfcommands`

#### BigWig Environment: `ribo_bigwig`
Required for strand-specific bigWig file generation:
```bash
conda env create -f ribo_bigwig_environment.yaml
```

This separate environment contains:
- `python=3.7`
- `deeptools>=3.5.0` (for bamCoverage)
- `samtools>=1.4` (for BAM handling)
- `numpy` and `pysam` (dependencies)

The pipeline automatically switches to `ribo_bigwig` when generating bigWig files. You only need to run the pipeline from the `ribo_genome` environment.

### HISAT2 Genome Index
Pre-built GRCh38 index files with Gencode v48 annotations:
- Generated using `hisat2_extract_exons.py` and `hisat2_extract_splice_sites.py`
- Stored under `./rf_sample_data/genome/`

---

## YAML Configuration

Your pipeline behavior is controlled through a `.yaml` configuration file. Below are the key sections you need to configure:

### Core Pipeline Settings

```yaml
# Disable ribo creation for genome alignment mode
do_ribo_creation: false

# Set deduplication method for Ribo-seq data
dedup_method: "umi_tools"  # Options: "none", "position", "umi_tools"
```

### Output Directory Configuration

```yaml
output:
  individual_lane_directory: 'individual'
  merged_lane_directory: 'merged'
  intermediates:
    base: 'intermediates_umi'
    clip: 'clip'
    log: 'log'
    transcriptome_alignment: 'transcriptome_alignment'
    filter: 'filter'
    genome_alignment: 'genome_alignment'
    bam_to_bed: 'bam_to_bed'
    quality_filter: 'quality_filter'
    deduplication: 'deduplication'
    bigwigs: 'bigwigs'  # BigWig files directory
  output:
    base: 'output_umi'
    log: 'log'
    fastqc: 'fastqc'
    ribo: 'ribo'
```

### P-site Offset Correction for Ribo-Seq Reads

Works with `dedup_method: "umi_tools"`, `"position"`, or `"none"`.

```yaml
psite_offset:
  offset_file: "./outputHuman_OFFSETS.csv"
  # Map each sample and lane to its source experiment ID from the offset file
  sample_matching:
    experiment_1:
      lanes:
        1: "GSM4114767"
        2: "GSM4114768"
    # Add more sample mappings as needed
```
### RNA-seq Configuration

```yaml
rnaseq:
  dedup_method: "position"  # Options: "none", "position"
  library_strandedness: "reverse"  # Options: "reverse" (default), "forward"
  hisat2_arguments: "--no-unal -k 1 --no-softclip"  # Genome alignment arguments

  # Optional read trimming settings
  genome_read_trim:
    enabled: true
    length: 35  # Trim reads to this length after rRNA filtering, before genome alignment
```

Note: `library_strandedness` controls how RNA-seq bigWig files are generated. Use `"reverse"` for most standard RNA-seq libraries (default), or `"forward"` if your library preparation uses forward-stranded protocols. `genome_read_trim` trims RNA-seq reads to a specified length (using cutadapt) after rRNA filtering and before genome alignment.

### Sample Input Configuration

```yaml
input:
  fastq_base: ""  # Optional base path for FASTQ files
  fastq:
    experiment_1:  # Sample/experiment name
      - "path/to/experiment_1_lane1.fastq.gz"
      - "path/to/experiment_1_lane2.fastq.gz"

rnaseq:
  fastq_base: ""  # Optional base path for RNA-seq FASTQ files
  fastq:
    experiment_1:  # Sample/experiment name (can be same or different from ribo-seq)
      - "path/to/experiment_1_rnaseq_lane1.fastq.gz"
      - "path/to/experiment_1_rnaseq_lane2.fastq.gz"
```

Create your configuration file based on these examples and adjust paths and parameters according to your dataset.

---
## P-site Offset Correction (Ribo-seq)

The pipeline performs **P-site offset correction** on ribo-seq data to extract the exact ribosome P-site position from each read.

### What's Generated

For each ribo-seq sample:
- P-site corrected BAM: `{sample}.psite.bam` - Each read converted to 1bp at P-site position
- P-site corrected BED: `{sample}.psite.bed` - BED format of P-site positions

### How It Works

1. Reads CSV file containing experiment-specific P-site offsets per read length
2. Maps each sample to an experiment ID via YAML configuration
3. Applies offsets based on read length:
   - Forward strand: Offset from 5' end
   - Reverse strand: Offset from 3' end
4. Outputs 1bp reads at exact P-site positions

### Requirements

- Deduplication method: Works with `dedup_method: "umi_tools"`, `"position"`, or `"none"`
- Sample mapping: Each sample and lane must be mapped to an experiment ID in YAML (using `sample_matching` with `lanes`)

---
## Strand-Specific BigWig File Generation

### What's Generated

Ribo-seq (from P-site corrected BAMs):
- Plus/Forward strand: `{sample}.psite.plus.bigWig`
- Minus/Reverse strand: `{sample}.psite.minus.bigWig`

RNA-seq (from deduplicated BAMs when `dedup_method: "position"`):
- Plus/Forward strand: `{experiment}.rnaseq.dedup.plus.bigWig`
- Minus/Reverse strand: `{experiment}.rnaseq.dedup.minus.bigWig`

RNA-seq (from quality-filtered BAMs when `dedup_method: "none"`):
- Plus/Forward strand: `{sample}.rnaseq.nodedup.plus.bigWig`
- Minus/Reverse strand: `{sample}.rnaseq.nodedup.minus.bigWig`

For RNA-seq with deduplication, samples are grouped by experiment name (removing the last dot and suffix, e.g., `experiment_1.1` → `experiment_1`) before merging for bigWig generation. The `library_strandedness` parameter controls how genomic strands are mapped to RNA strands in the bigWig files.

### Output Directory Structure

```
intermediates_umi/
├── deduplication/
│   ├── merged/
│   │   ├── {sample}.dedup.bam                    # Ribo-seq: UMI-deduplicated
│   │   ├── {sample}.dedup.bed
│   │   ├── {sample}.psite.bam                    # Ribo-seq: P-site corrected
│   │   └── {sample}.psite.bed
│   │
│   └── bigwigs/merged/                           # Ribo-seq bigWigs
│       ├── {sample}.psite.plus.bigWig
│       └── {sample}.psite.minus.bigWig
│
└── rnaseq/
    └── deduplication/
        ├── merged/
        │   ├── {sample}.rnaseq_genome.post_dedup.bed
        │   ├── {sample}.rnaseq_genome.post_dedup.bam
        │   └── {sample}.rnaseq_genome.post_dedup.bam.bai
        │
        └── bigwigs/merged/                       # RNA-seq bigWigs
            ├── {experiment}.rnaseq.dedup.plus.bigWig    # With dedup (position)
            ├── {experiment}.rnaseq.dedup.minus.bigWig
            ├── {sample}.rnaseq.nodedup.plus.bigWig      # Without dedup (none)
            └── {sample}.rnaseq.nodedup.minus.bigWig
```

---

## Directory Structure Overview

### Ribo-seq Pipeline
```
intermediates_umi/
├── clip/                              # Adapter clipping
├── umi_tools/merged/                  # UMI extraction
├── filter/                            # rRNA filtering
├── genome_alignment/
│   ├── individual/                    # Per-lane alignments
│   └── merged/                        # Merged alignments
├── quality_filter/                    # Quality filtered alignments
├── bam_to_bed/individual/             # BED conversion
├── deduplication/merged/              # All deduplication outputs
│   ├── *.dedup.bam                    # UMI-deduplicated
│   ├── *.psite.bam                    # P-site corrected
│   └── *.bed files
└── deduplication/bigwigs/merged/      # Coverage tracks
    └── *.psite.{plus,minus}.bigWig
```

### RNA-seq Pipeline
```
intermediates_umi/rnaseq/
├── clip/                              # Adapter clipping
├── filter/                            # rRNA filtering
├── genome_alignment/
│   ├── individual/                    # Per-lane alignments
│   └── merged/                        # Merged alignments
├── quality_filter/                    # Quality filtered alignments
├── bam_to_bed/
│   ├── individual/                    # Individual BED files
│   └── merged/                        # Merged BED files
├── deduplication/merged/              # All deduplication outputs
│   ├── *.post_dedup.bed               # Position-deduplicated
│   └── *.post_dedup.bam               # For BigWig generation
└── deduplication/bigwigs/merged/      # Coverage tracks
    ├── *.rnaseq.dedup.{plus,minus}.bigWig    # With dedup (position)
    └── *.rnaseq.nodedup.{plus,minus}.bigWig  # Without dedup (none)
```

---

## Notes

P-site corrected BAM files contain 1bp reads positioned at the ribosome's P-site.

The `library_strandedness` parameter controls how RNA-seq reads are assigned to strands in bigWig files. Use `"reverse"` (default) for most standard RNA-seq libraries where the first read maps to the opposite strand of the gene. Use `"forward"` if your library preparation protocol produces forward-stranded libraries. This parameter only affects RNA-seq bigWig generation, not ribo-seq. For forward-stranded libraries (`++`, `--`), a read mapped to the '+' strand indicates the parental gene is on the '+' strand, and a read mapped to the '-' strand indicates the parental gene is on the '-' strand. For reverse-stranded libraries (`+-`, `-+`), a read mapped to the '-' strand indicates the parental gene is on the '+' strand. You can use `infer_experiment.py` from the RSeQC package to determine the strandedness of your RNA-seq library by analyzing the orientation of reads relative to gene features.

---

## Version Information

### Current Version: 2.0.0

### Recent Changes
- Added strand-specific BigWig generation
- P-site offset correction
- Library strandedness parameter
- RNA-Seq Trimming
- Statistics compilation with genome alignment support
- Genome alignment support

### Compatibility
- Nextflow: 19.04.1
- Conda: 4.12+
- Python: 3.7+
- Java: OpenJDK 8+

---

## Additional Documentation

For more detailed information, see:

- **[Code Audit Report](docs/CODE_AUDIT.md)**: Comprehensive audit of all processes, channels, functions, commands, and scripts
- **[Changes from Reference Version](docs/RiboFlow_Genome_Changes.md)**: Detailed comparison with the original RiboFlow pipeline
- **[Scripts Documentation](scripts/README.md)**: Complete reference for all Python utility scripts
- **[Commands Documentation](modified_rfcommands/README.md)**: Complete reference for all `rfc` command-line tools

## Code Structure

The pipeline consists of several key components:

### Main Pipeline
- **`RiboFlow.groovy`**: Main Nextflow pipeline script (3,253 lines)
  - 42 processes for Ribo-seq and RNA-seq processing
  - Comprehensive channel architecture for data flow
  - Conditional execution based on configuration

### Command-Line Tools
- **`modified_rfcommands/`**: Custom `rfc` commands
  - 7 actively used commands
  - HISAT2 log processing
  - Statistics compilation and merging
  - Deduplication utilities

### Utility Scripts
- **`scripts/`**: Python utility scripts
  - 6 scripts for specialized processing
  - P-site offset correction
  - Read extraction
  - Statistics updates

### Configuration
- **YAML files**: Pipeline configuration
  - `project_umi.yaml`: UMI deduplication example
  - `project_position.yaml`: Position deduplication example
  - `project_nodedup.yaml`: No deduplication example

