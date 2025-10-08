# RiboFlow_Genome
---

## Quick Start

```bash
# 1. Create conda environments (one-time setup)
conda env create -f environment.yaml
conda env create -f ribo_bigwig_environment.yaml

# 2. Activate main environment
conda activate ribo_genome

# 3. Run the pipeline
nextflow run RiboFlow.groovy -params-file project_umi.yaml # -c configs/local.config (optional)
```

---
## Installation & Dependencies

### **Conda Environments**

#### **1. Main Environment: `ribo_genome`**
Created from `environment.yaml` - contains all core pipeline tools:
```bash
conda env create -f environment.yaml
```

#### **2. BigWig Environment: `ribo_bigwig`**
Required for strand-specific bigWig file generation:
```bash
conda env create -f ribo_bigwig_environment.yaml
```

This separate environment contains:
- `deeptools` (for bamCoverage)
- `samtools` (for BAM handling)
- `numpy` and `pysam` (dependencies)

**Note**: The pipeline automatically switches to `ribo_bigwig` when generating bigWig files. You only need to run the pipeline from the `ribo_genome` environment.

### **HISAT2 Genome Index**
Pre-built GRCh38 index files with Gencode v48 annotations:
- Generated using `hisat2_extract_exons.py` and `hisat2_extract_splice_sites.py`
- Stored under `./rf_sample_data/genome/`

---

## YAML Configuration Changes

### **Core Pipeline Settings**
```yaml

# ribo_creation should be set to false, as genome alignments aren't supported.
do_ribo_creation: false

# Independent deduplication methods
dedup_method: "umi_tools"           # Ribo-seq: "none", "position", "umi_tools"


# Output folder settings
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
    bigwigs: 'bigwigs'                # BigWig files directory
  output:
    base: 'output_umi'
    log: 'log'
    fastqc: 'fastqc'
    ribo: 'ribo'
```

### **P-site Offset Correction (Optional)**
```yaml
psite_offset:
  offset_file: "./outputHuman_OFFSETS.csv"
  # Map each ribo-seq sample to its source experiment ID from the offset file
  experiment_mapping:
    1cell-2: "GSM3907597"
    1cell-4: "GSM3984249"
```

### **RNA-seq Genome Alignment**
```yaml
rnaseq:
  dedup_method: "position"          # RNA-seq: "none", "position"
  hisat2_arguments: "--no-unal -k 1 -L 15" # genome arguments

  # Optional read trimming
  transcriptome_read_trim:
    enabled: true
    length: 35

  genome_read_trim:
    enabled: true
    length: 35
```

---
## P-site Offset Correction (Ribo-seq)

The pipeline performs **P-site offset correction** on ribo-seq data to extract the exact ribosome P-site position from each read.

### **What's Generated**

For each ribo-seq sample:
- **P-site corrected BAM**: `{sample}.psite.bam` - Each read converted to 1bp at P-site position
- **P-site corrected BED**: `{sample}.psite.bed` - BED format of P-site positions

### **How It Works**

1. Reads CSV file containing experiment-specific P-site offsets per read length
2. Maps each sample to an experiment ID via YAML configuration
3. Applies offsets based on read length:
   - **Forward strand**: Offset from 5' end
   - **Reverse strand**: Offset from 3' end
4. Outputs 1bp reads at exact P-site positions

### **Requirements**

- **Deduplication method**: Only works with `dedup_method: "umi_tools"`
- **Experiment mapping**: Each sample must be mapped to an experiment ID in YAML

### **Output Directory**
```
intermediates_umi/deduplication/merged/
├── {sample}.dedup.bam              # UMI-deduplicated BAM
├── {sample}.dedup.bam.bai
├── {sample}.psite.bam              # P-site corrected BAM (1bp reads)
├── {sample}.psite.bam.bai
└── {sample}.psite.bed              # P-site positions in BED format
```

---
## Strand-Specific BigWig File Generation
### **What's Generated**

**Ribo-seq** (from P-site corrected BAMs):
- **Plus/Forward strand**: `{sample}.psite.plus.bigWig`
- **Minus/Reverse strand**: `{sample}.psite.minus.bigWig`

**RNA-seq** (from deduplicated BAMs):
- **Plus/Forward strand**: `{sample}.rnaseq.dedup.plus.bigWig`
- **Minus/Reverse strand**: `{sample}.rnaseq.dedup.minus.bigWig`

### **Output Directory Structure**

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
            ├── {sample}.rnaseq.dedup.plus.bigWig
            └── {sample}.rnaseq.dedup.minus.bigWig
```

### **Technical Details**

- **Ribo-seq**: Generated from P-site corrected BAMs (requires UMI deduplication + P-site correction)
- **RNA-seq**: Generated from position-deduplicated BAMs (requires `rnaseq.dedup_method == 'position'`)
- **Resolution**: 1bp bins (`--binSize 1`)
- **Strand filtering**: Uses `bamCoverage --filterRNAstrand` for strand separation

### **Processes Added**

**Ribo-seq**:
1. **`apply_psite_correction`**: Converts deduplicated reads to 1bp P-site positions
2. **`genome_create_strand_specific_bigwigs`**: Generates strand-specific bigWigs from P-site BAMs

**RNA-seq**:
1. **`rnaseq_genome_deduplicate`**: Position-based deduplication on merged BED
2. **`rnaseq_dedup_bed_to_bam`**: Converts deduplicated BED to BAM
3. **`rnaseq_create_strand_specific_bigwigs`**: Generates strand-specific bigWigs

---

## Directory Structure Overview

### **Ribo-seq Pipeline**
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

### **RNA-seq Pipeline**
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
    └── *.rnaseq.dedup.{plus,minus}.bigWig
```
---
