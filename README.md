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
* [Frequently Asked Questions](https://github.com/ribosomeprofiling/riboflow/blob/master/FAQ.md)  
* [Release Notes](https://github.com/ribosomeprofiling/riboflow/blob/master/CHANGELOG.md)  

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
# NOTE: The sample data does not contain the HISAT2 indexes for genome alignment as they are too large.
# Please contact daniel_nguyen@utexas.edu if you need the indexes.
# Additionally, rf_sample_data does not yet contain @rf_sample_data/annotation/hg38_GENCODE_V47.bed.gz
# needed for automatic library strandedness detection.

nextflow RiboFlow.groovy -params-file project_nodedup.yaml -profile docker_local
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
# NOTE: The sample data does not contain the HISAT2 indexes for genome alignment as they are too large.
# Please contact daniel_nguyen@utexas.edu if you need the indexes.
# Additionally, rf_sample_data does not yet contain @rf_sample_data/annotation/hg38_GENCODE_V47.bed.gz
# needed for automatic library strandedness detection.

# Finally run RiboFlow
nextflow RiboFlow.groovy -params-file project_nodedup.yaml

```

## Output

Pipeline run may take several minutes.
When finished, the resulting files are in the `./output` folder.

Mapping statistics are compiled in a csv file called `stats.csv`

```
ls output/stats/stats.csv
```

### BigWig Files

BigWig files are generated for both Ribo-Seq and RNA-Seq data.

*   **Ribo-Seq BigWigs**: Created in `intermediates/alignment_ribo/bigwigs`.
    *   **Note:** These represent P-sites (requires P-site correction to be enabled).
*   **RNA-Seq BigWigs**: Automatically generated in `intermediates/rnaseq/alignment_ribo/bigwigs/merged/`.
    *   **Note:** These represent the coverage of the length of aligned reads (whether trimmed or not).

BigWig generation does not require P-site correction for RNA-Seq, but Ribo-Seq bigwigs are based on P-sites.


## RiboFlow on Your Data

For running RiboFlow on actual data, files must be organized and a parameters file must be prepared.
You can examine the sample run above to see an example.

1. Organize your data. The following files are required for RiboFlow
   * **Ribosome profiling sequencing data:** in gzipped fastq files
   * **Transcriptome Reference:** Bowtie2 index files
   * **Filter Reference:** Bowtie2 index files (typically for rRNA sequences)
   * **Annotation:** A bed file defining CDS, UTR5 and UTR3 regions.
   * **Transcript Lengths:** A two column tsv file containing transcript lengths

2. Prepare a custom `project.yaml` file.
You can use the sample file `project_nodedup.yaml`, provided in this repository,
as template.

3. In `project.yaml`, provide RiboFlow parameters such as `clip_arguments`, alignment arguments etc.
You can simply modify the arguments in the sample file `project_nodedup.yaml` in this repository.

4. You can adjust the hardware and computing environment settings in Nextflow configuration file(s).
For Docker option, see `configs/docker_local.config`. If you are not using Docker,
see `configs/local.config`.

5. RNA-Seq data is optional for RiboFlow. So, if you do NOT have RNA-Seq data, in the project file, set

`do_rnaseq: false`

If you have RNA-Seq data to be paired with ribosome profiling data, see the [Advanced Features](#advanced-features) below.


6. Metadata is optional for RiboFlow. If you do NOT have metadata, in the project file, set

`do_metadata: false`

If you have metadata, see [Advanced Features](#advanced-features) below.

7. Run RiboFlow using the new parameters file `project.yaml`.

Using Docker:

`nextflow RiboFlow.groovy -params-file project_nodedup.yaml -profile docker_local`

Without Docker:

`nextflow RiboFlow.groovy -params-file project_nodedup.yaml`

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
# NOTE: The sample data does not contain the HISAT2 indexes for genome alignment as they are too large.
# Please contact daniel_nguyen@utexas.edu if you need the indexes.
# Additionally, rf_sample_data does not yet contain @rf_sample_data/annotation/hg38_GENCODE_V47.bed.gz
# needed for automatic library strandedness detection.

# Finally run RiboFlow
nextflow RiboFlow.groovy -params-file project_umi.yaml
```

 ### UMI support for RNA-Seq
 
 In the current version, UMIs are supported for ribosome profiling data only. So RNA-Seq libraries can either be used without deduplication or the reads can be collapsed based on position.

## A Note on References

**RiboFlow_genome** is designed to work with **genomic references**. This version is **NOT** configured for transcriptome alignment.

For transcriptome alignment, please use the original [RiboFlow](https://github.com/ribosomeprofiling/riboflow).

The users need to provide a genomic reference and annotation to run this software.

## Advanced Features

### P-site Offset Correction
You can enable P-site offset correction by configuring `psite_offset` in your project YAML. This requires a CSV file with experiment-specific P-site offsets per read length.
See `project_nodedup.yaml` for an example configuration.

### Library Strandedness
You can specify the library strandedness in the configuration using `library_strandedness`.
Options are `"forward"`, `"reverse"`, or `"auto"`. If set to `"auto"`, the pipeline will attempt to detect strandedness automatically (requires `hg38_GENCODE_V47.bed.gz` annotation).
See `project_nodedup.yaml` for configuration details.

### Genome Read Trim
You can optionally trim reads after filtering but before alignment using the `genome_read_trim` setting.
See `project_nodedup.yaml` for details.

### RNA-Seq Data

If you have RNA-Seq data that you want to pair with ribosome profiling experiments,
provide the paths of the RNA-Seq (gzipped) fastq files  in the configuration file in
_input -> metadata_. See the file `project_nodedup.yaml` in this repository for an example.
Note that the names in defining RNA-Seq files must match the names in definig ribosome profiling data.
Also turn set the do_rnaseq flag to true, in the project file:

`do_rnaseq: true`

Transcript abundance data will be processed.

### Metadata

If you have metadata files for the ribosome profiling experiments,
provide the paths of the metadata files (in yaml format) in the configuration file in
_input -> metadata_. See the file `project_nodedup.yaml` in this repository for an example.
Note that the names in defining metadata files must match the names in definig ribosome profiling data.
Also turn set the metadata flag to true, in the project file:

`do_metadata: true`

Metadata will be processed.

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
