[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3376949.svg)](https://doi.org/10.5281/zenodo.3376949)

![RiboFlow](/docs/figures/riboflow_logo.jpg "RiboFlow Logo")

# RiboFlow_genome #

RiboFlow_genome is a forked version of the original [RiboFlow](https://github.com/ribosomeprofiling/riboflow) pipeline, customized for genomic alignment.
RiboFlow_genome belongs to a [software ecosystem](https://ribosomeprofiling.github.io/) desgined to work with ribosome profiling data.


![Overview](/docs/figures/ecosystem_overview.jpg "Ribo Ecosystem Overview")

## Contents

* [Installation](#installation)
* [Test Run](#test-run)
* [Output](#output)
* [RiboFlow on Your Data](#riboflow-on-your-data)
* [UMIs](#working-with-unique-molecular-identifiers)
* [A Note on References](#a-note-on-references)
* [Advanced Features](#advanced-features)
* [FAQ](FAQ.md)
* [Changelog](CHANGELOG.md)

## Installation

### Requirements

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://docs.docker.com/install/) (Optional)
* [Conda](https://conda.io/en/latest/miniconda.html) (Optional)

First, follow the instructions in [Nextflow website](https://www.nextflow.io/) and install Nextflow.

### Docker Option

Install [Docker](https://docs.docker.com/install/).
Here is a [tutorial for Ubuntu.](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04)

All remaining dependencies come in the Docker image [danielnguyener/riboflow](https://hub.docker.com/r/danielnguyener/riboflow).
This image is automatically pulled by RiboFlow when run with Docker (see test runs below).

### Conda Option

This option has been tested on Linux systems only.

Install  [Conda](https://conda.io/en/latest/miniconda.html).

All other dependencies can be installed using the environment file,
environment.yaml, in this repository.
```
git clone https://github.com/DanielNguyener/riboflow_genome.git
conda env create -f riboflow_genome/environment.yaml
```

The above command will create a conda environment called _ribo_genome_
and install dependencies in it.
To start using RiboFlow, you need to activate the _ribo_genome_ environment.

`conda activate ribo_genome`

## Test Run

For fresh installations, before running RiboFlow on actual data,
it is recommended to do a test run.


### Run Using Docker

```
# Clone this repository in a new folder and change your working directory to the RiboFlow folder.
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/DanielNguyener/riboflow_genome.git
cd riboflow_genome

# Obtain a copy of the sample data in the working directory.
git clone https://github.com/ribosomeprofiling/rf_sample_data.git
# NOTE: The sample data does NOT contain the STAR genome index (too large to ship).
# Build one yourself with `STAR --runMode genomeGenerate` against your reference of
# choice, and point input.reference.genome at the resulting directory in your params file.

nextflow RiboFlow.groovy -params-file example.yaml -profile docker_local
```

Note that we provided the argument `-profile docker_local` to Nextflow to indicate that RiboFlow will be run via Docker containers. In other words, the steps of RiboFlow will be executed inside Docker containers by Nextflow. 
Hence, no locally installed software (other than Java and Nextflow) is needed by RiboFlow.  


### Run Using Conda Environment

In Conda option, the steps of RiboFlow are run locally. So, we need to install the dependencies first. This can easily be done via conda. The default profile directs RiboFlow to run locally, so we can simply skip the `-profile` argument. Also note that the conda environment has to be activated before running RiboFlow. 

Before running the commands below, make sure that you have created the conda environment, called _ribo_genome_,
using the instructions above. 

```
# List the environments to make sure that ribo_genome environment exists
conda env list

# Activate the ribo_genome environment
conda activate ribo_genome

# Get RiboFlow repository
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/DanielNguyener/riboflow_genome.git
cd riboflow_genome

# Obtain a copy of the sample data in the working directory.
git clone https://github.com/ribosomeprofiling/rf_sample_data.git
# NOTE: The sample data does NOT contain the STAR genome index (too large to ship).
# Build one yourself with `STAR --runMode genomeGenerate` and point
# input.reference.genome at the resulting directory in your params file.

# Finally run RiboFlow
nextflow RiboFlow.groovy -params-file example.yaml

```

## Output

A successful run produces this layout:

| Path | Contents |
|---|---|
| `output/stats/stats.csv` | Per-sample alignment summary (one row per sample) |
| `output/stats/individual_stats.csv` | Per-lane alignment summary (one row per lane) |
| `output/bigwigs/ribo/*.ribo.{plus,minus}.bigWig` | Strand-specific ribo-seq bigWigs |
| `output/bigwigs/rnaseq/*.rnaseq.bigWig` | RNA-seq coverage bigWigs (unstranded) |
| `output/alignments/ribo/individual/*.bam, *.bed` | Per-lane deduplicated ribo-seq alignments |
| `output/alignments/ribo/merged/*.bam, *.bed` | Merged-sample deduplicated ribo-seq alignments |
| `output/alignments/rnaseq/{individual,merged}/` | RNA-seq deduplicated alignments |
| `output/rnaseq/stats/` | RNA-seq alignment stats CSVs |
| `output/fastqc/` | FastQC reports (if `do_fastqc: true`) |
| `intermediates/` | Cached working files (raw STAR BAMs, qpass BAMs, pre-dedup BEDs). Safe to delete; will be regenerated on re-run. |

Ribo-seq bigwigs cover read 5' ends on the genome (no P-site correction; this is a pure
genome-alignment pipeline). RNA-seq bigwigs are unstranded coverage.

### Compute requirements

Right-sizing depends on the genome and dataset, but a useful baseline:

- **STAR genome alignment**: ~30 GB RAM for the human genome index, ~10–15 minutes per
  ribo-seq lane on 16 cores.
- **Bowtie2 rRNA filter**: ~2 GB RAM, single-thread-bound for short reads.
- **umicollapse dedup**: ~32 GB JVM heap on a CCLE-scale BAM. Memory-bound.
- **Disk**: budget ~3-5x the input FASTQ size for `intermediates/`.

See `configs/local.config` for a 128-core / 240 GB workstation profile and
`configs/example_slurm_cluster.config` for a SLURM-style starting point.


## RiboFlow on Your Data

For running RiboFlow on actual data, files must be organized and a parameters file must be prepared.
You can examine the sample run above to see an example.

1. Organize your data. The following are required:
   * **Ribosome profiling sequencing data:** gzipped FASTQ files (single-end).
   * **(Optional) RNA-seq sequencing data:** gzipped FASTQ, single- or paired-end.
   * **Filter Reference:** a Bowtie2 index for rRNA / contaminant filtering. The
     upstream [references_for_riboflow](https://github.com/ribosomeprofiling/references_for_riboflow)
     repository has working filter indices for several organisms.
   * **Genome Reference:** a STAR genome index directory (built with the same STAR
     version pinned in `environment.yaml`).

2. Prepare a custom `project.yaml` file. Copy `example.yaml` as a template.

3. In `project.yaml`, provide parameters such as `clip_arguments`, STAR alignment
   arguments, dedup method, etc. See the comments in `example.yaml` for the full
   set of options.

4. You can adjust the hardware and computing environment settings in Nextflow configuration file(s).
For Docker option, see `configs/docker_local.config`. If you are not using Docker,
see `configs/local.config`.

5. RNA-Seq data is optional for RiboFlow. So, if you do NOT have RNA-Seq data, in the project file, set

`do_rnaseq: false`

If you have RNA-Seq data to be paired with ribosome profiling data, see the [Advanced Features](#advanced-features) below.


6. Run RiboFlow using the new parameters file `project.yaml`.

Using Docker:

`nextflow RiboFlow.groovy -params-file project.yaml -profile docker_local`

Without Docker:

`nextflow RiboFlow.groovy -params-file project.yaml`

## Working with Unique Molecular Identifiers
Unique Molecular Identifiers (UMIs) can be ligated to either side of the molecules and 
they allow labeling molecules uniqely. This way UMIs can be used to deduplicate mapped reads
for more accurate quantification.

If there are UMIs in your ribosome profiling data, Riboflow can trim them and deduplicate reads based on UMIs. 

RiboFlow extracts UMIs and stores them in the Fastq headers and uses the UMIs in deduplication 
(instead of position based read collapsing). For this purpose RiboFlow uses 
[umi_tools](https://github.com/CGATOxford/UMI-tools).


### Project File

Here we explain the related parts of the project file to be able to use UMIs feature of Riboflow.

Also, we provide a working example of project file in this repository: *project_umi.yaml*.

The following parameter must be set:
```
dedup_method: "umi_tools"
```

Also, users must set the following two parameters: `umi_tools_extract_arguments` and `umi_tools_dedup_arguments`.

For example: 
```
umi_tools_extract_arguments: "-p \"^(?P<umi_1>.{12})(?P<discard_1>.{4}).+$\" --extract-method=regex"
umi_tools_dedup_arguments:   "--read-length"
```

The above example takes the first 12 nucleotides from the 5' end, discards the 4 nucleotides downstream and writes the 12 nt UMI sequence to the header.
The second parameter tells umi_tools to use read lengths IN ADDITION to UMI sequencing in collapsing reads. Note that these two arguments are directly provided to umi_tools. So users are encouraged to familirize themselves with [umi_tools](https://umi-tools.readthedocs.io/en/latest/).

### Test Run with UMIs

We provide a mini dataset, with two samples, to try Riboflow with sequencing reads having UMIs.
In this sample dataset, the first 12 nucleotides on the 5' end of the reads are UMIs.
Four nucleotides following the UMIs need to be discarded.
On the 3' end of the reads, there are adapters having the sequence `AAAAAAAAAACAAAAAAAAAA`.
The parameters of this sample run are provided in the file *project_umi.yaml*.
Below are the steps to process this data.


```
# List the environments to make sure that ribo_genome environment exists
conda env list

# Activate the ribo_genome environment
conda activate ribo_genome

# Get RiboFlow repository
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/DanielNguyener/riboflow_genome.git
cd riboflow_genome

# Obtain a copy of the sample data in the working directory.
git clone https://github.com/ribosomeprofiling/rf_sample_data.git
# NOTE: The sample data does NOT contain the STAR genome index (too large to ship).
# Build one yourself with `STAR --runMode genomeGenerate` and point
# input.reference.genome at the resulting directory in your params file.

# Finally run RiboFlow
nextflow RiboFlow.groovy -params-file example.yaml
```

### UMI support for RNA-seq

UMIs are supported for ribo-seq data via `dedup_method: "umicollapse"`. RNA-seq lanes
either skip dedup or use coordinate-based dedup (`rnaseq.dedup_method: "position"`).

## A Note on References

**riboflow_genome** is designed to work with **genomic references** (STAR index +
Bowtie2 rRNA filter). This version is **NOT** configured for transcriptome alignment.

For transcriptome alignment, please use the original
[RiboFlow](https://github.com/ribosomeprofiling/riboflow).

## Advanced Features

### Pairing Ribo-seq with RNA-seq

Set `do_rnaseq: true` in your params file and provide RNA-seq FASTQs under the
`rnaseq.fastq.<sample>` keys. Sample names must match the ribo-seq sample names.
RNA-seq lanes can be either single-end (string) or paired-end (`[R1, R2]` list).
See `example.yaml` for the full schema.

### Tuning quality-pass filters

`mapping_quality_cutoff` (ribo) and `rnaseq.mapping_quality_cutoff` set the MAPQ
threshold for the qpass step. `ribo_filter_flags` (default `2052`) and
`rnaseq.filter_flags` (default `2308`) are the SAM flag masks passed to
`samtools view -F`. The ribo default keeps secondary alignments so multi-mappers
contribute to bigwigs; the RNA-seq default drops them to enforce unique-only
output.

## Citing

[RiboFlow, RiboR and RiboPy: an ecosystem for analyzing ribosome profiling data at read length resolution, H. Ozadam, M. Geng, C. Cenik
Bioinformatics 36 (9), 2929-2931](https://academic.oup.com/bioinformatics/article/36/9/2929/5701654)

```bibtex
@article{ozadam2020riboflow,
  title={RiboFlow, RiboR and RiboPy: an ecosystem for analyzing ribosome profiling data at read length resolution},
  author={Ozadam, Hakan and Geng, Michael and Cenik, Can},
  journal={Bioinformatics},
  volume={36},
  number={9},
  pages={2929--2931},
  year={2020},
  publisher={Oxford University Press}
}
```

## [Frequently Asked Questions](https://github.com/ribosomeprofiling/riboflow/blob/master/FAQ.md)  

  
## [Release Notes](https://github.com/ribosomeprofiling/riboflow/blob/master/CHANGELOG.md)  
