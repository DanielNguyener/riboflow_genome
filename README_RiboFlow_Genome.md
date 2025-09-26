# RiboFlow_Genome
---
## Installation & Dependencies

### **Modified rfcommands Package**
The enhanced pipeline requires updated rfcommands with new HISAT2 log processing capabilities:

```bash
# Remove old installation
pip uninstall rfcommands-riboflow -y
rm -rf ~/.local/lib/python*/site-packages/rfcommands*
rm -rf ./modified_rfcommands/rfcommands_riboflow.egg-info/

# Clean install
cd ./modified_rfcommands/
pip install -e .

# Verify installation
rfc hisat2-log-to-csv --help
rfc merge-hisat2-logs --help
```

### **HISAT2 Genome Index**
Pre-built GRCh38 index files with Gencode v48 annotations:
- Generated using `hisat2_extract_exons.py` and `hisat2_extract_splice_sites.py`
- Stored under `./rf_sample_data/genome/`

---

## YAML Configuration Changes

### **Core Pipeline Settings**
```yaml
# Only supports transcriptome alignments.
do_ribo_creation: true

# Independent deduplication methods
dedup_method: "umi_tools"           # Ribo-seq: "none", "position", "umi_tools"

# Comment out to skip transcriptome or genome alignment
input:
  reference:
    transcriptome: ./rf_sample_data/transcriptome/appris_human_24_01_2019_selected*
    genome: ./rf_sample_data/genome/hisat2_spliced_index*
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
