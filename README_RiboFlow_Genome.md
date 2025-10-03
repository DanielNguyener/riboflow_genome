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
nextflow run RiboFlow.groovy -params-file project_umi.yaml -c configs/local.config
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
    bigwigs: 'bigwigs'                # NEW: BigWig files directory
  output:
    base: 'output_umi'
    log: 'log'
    fastqc: 'fastqc'
    ribo: 'ribo'
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
## Strand-Specific BigWig File Generation

The pipeline generates **strand-specific bigWig files** from deduplicated RNA-seq data for genome browser visualization.

### **What's Generated**

For each RNA-seq sample, two bigWig files are created per strand:
- **Plus/Forward strand**: `{sample}.rnaseq.dedup.plus.bigWig`
- **Minus/Reverse strand**: `{sample}.rnaseq.dedup.minus.bigWig`

BigWigs are generated at **two levels**:
1. **Individual lanes**: Each sample/lane gets its own pair of bigWig files
2. **Merged experiments**: Combined lanes per experiment get merged bigWig files

### **Output Directory Structure**

```
intermediates_umi/rnaseq/deduplication/
├── merged/
│   ├── {sample}.rnaseq_genome.post_dedup.bed     
│   ├── {sample}.rnaseq_genome.post_dedup.bam     
│   └── {sample}.rnaseq_genome.post_dedup.bam.bai 
│
├── individual/
│   ├── {sample}.{index}.rnaseq_genome.post_dedup.bed
│   ├── {sample}.{index}.rnaseq_genome.post_dedup.bam      
│   └── {sample}.{index}.rnaseq_genome.post_dedup.bam.bai  
│
└── bigwigs/                                        
    ├── individual/                                 # Per-sample bigWigs
    │   ├── {sample}.{index}.rnaseq.dedup.plus.bigWig
    │   └── {sample}.{index}.rnaseq.dedup.minus.bigWig
    │
    ├── {sample}.rnaseq.dedup.plus.bigWig          # Per-experiment bigWigs
    └── {sample}.rnaseq.dedup.minus.bigWig
```

### **Technical Details**

- **Prerequisite**: Only runs when `rnaseq.dedup_method == 'position'`

### **New Processes Added**

1. **`rnaseq_individual_dedup_bed_to_bam`**: Converts individual deduplicated BED files to BAM
2. **`rnaseq_individual_create_strand_specific_bigwigs`**: Generates individual lane bigWigs
3. **`rnaseq_dedup_bed_to_bam`**: Converts merged deduplicated BED files to BAM
4. **`rnaseq_create_strand_specific_bigwigs`**: Generates merged experiment bigWigs

---
