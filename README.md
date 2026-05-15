[DOI](https://doi.org/10.5281/zenodo.3376949)

RiboFlow

# RiboFlow_genome

RiboFlow_genome is a forked version of the original [RiboFlow](https://github.com/ribosomeprofiling/riboflow) pipeline, customized for genomic alignment.
RiboFlow_genome belongs to a [software ecosystem](https://ribosomeprofiling.github.io/) desgined to work with ribosome profiling data.

Overview

## Contents

- [Installation](#installation)
  - [Four ways to run](#four-ways-to-run-local-or-hpc-conda-or-container)
  - [Docker Desktop (local)](#docker-desktop-local-workstations)
- [Test Run](#test-run)
- [Running on an HPC Cluster (Apptainer / Singularity)](#running-on-an-hpc-cluster-apptainer--singularity)
- [Output](#output)
- [Building the STAR genome index](#building-the-star-genome-index)
- [RiboFlow on Your Data](#riboflow-on-your-data)
- [UMIs](#working-with-unique-molecular-identifiers)
- [A Note on References](#a-note-on-references)
- [Advanced Features](#advanced-features)
- [FAQ](FAQ.md)
- [Changelog](CHANGELOG.md)

## Installation

### Requirements

- [Nextflow](https://www.nextflow.io/) **19.04.1** for running the DSL1 pipeline.
Modern Nextflow (≥22) only parses DSL1 in config-check mode; actual runs need
the 19.04.1 binary. Nextflow 19.04.1 requires Java 8 or 11.
- Choose **one** toolchain using the [four supported patterns](#four-ways-to-run-local-or-hpc-conda-or-container) (local vs HPC, conda vs container):
  - [Conda](https://conda.io/en/latest/miniconda.html) (Linux — local workstation **or** HPC node).
  - [Docker Desktop](https://www.docker.com/products/docker-desktop/)
  (`linux/amd64` interactive container — usual path on **macOS/Windows**, also works on Linux;
  [Docker Desktop](#docker-desktop-local-workstations)).
  - [Apptainer / Singularity](https://apptainer.org/) (HPC **container** shells;
  [HPC section](#running-on-an-hpc-cluster-apptainer--singularity)).

### Four ways to run (local or HPC, conda or container)

The pipeline invocation is always the same **inside** whichever environment has
STAR, Bowtie2, Java, Nextflow 19.x, … on `PATH`; only the *shell you start in*
changes:


| Pattern                         | Bring tools onto `PATH`                                                                                   | Notes                                                                                                                                                         |
| ------------------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Local • conda**               | Create/activate `**ribo_genome`** from `environment.yaml`.                                                | **Linux workstations only** — pinned Bioconda packages do **not** build on macOS. Step-by-step: [Conda option](#conda-option-linux-workstation-or-hpc-node).  |
| **Local • Docker**              | `docker pull/run --platform linux/amd64`, interactive `bash` inside the image, repo mounted into `/work`. | **macOS / Windows / Linux.** Image is **amd64**; Apple Silicon needs `--platform`. Details: [Docker Desktop](#docker-desktop-local-workstations).             |
| **HPC • conda** (e.g. TACC LS6) | Your Miniforge/Mamba `**ribo_genome`** on `$WORK`; `source …/conda.sh && conda activate ribo_genome`.     | Same `nextflow run` as everywhere else. Convenient on login/interactive nodes.                                                                                |
| **HPC • Apptainer**             | `apptainer pull docker://danielnguyener/riboflow:latest` then `**apptainer shell riboflow_latest.sif`**.  | **Recommended on Lustre** vs per-task `apptainer exec` (squashfuse / PATH issues). Details: [HPC section](#running-on-an-hpc-cluster-apptainer--singularity). |


**Nextflow command** (from the repo directory, tools already on `PATH`):

```bash
nextflow run RiboFlow.groovy -params-file project.yaml -c configs/local.config
```

Use `example.yaml` / `example_local.yaml` while learning. There is **no**
per-task Docker/Apptainer Nextflow profile in this repository — Docker and
Apptainer instructions mean **interactive container shells only**; the workflow
still uses profile `**standard`** (`configs/local.config`).

### Conda option (Linux workstation or HPC node)

```bash
git clone https://github.com/DanielNguyener/riboflow_genome.git
conda env create -f riboflow_genome/environment.yaml
conda activate ribo_genome
```

The environment is named `ribo_genome` and ships every binary the
pipeline shells out to (STAR, bowtie2, samtools, cutadapt, umicollapse,
deeptools, bedtools, java11, …). Activate it before invoking Nextflow.

The conda env is **Linux-only** — several pinned packages don't build
on macOS (use [Docker Desktop](#docker-desktop-local-workstations) there).
On **TACC or other Linux HPC**, conda is a first-class option alongside
**[Apptainer](#running-on-an-hpc-cluster-apptainer--singularity)** — see
[Four ways to run](#four-ways-to-run-local-or-hpc-conda-or-container).

### Docker Desktop (local workstations)

`[danielnguyener/riboflow:latest](https://hub.docker.com/r/danielnguyener/riboflow)`
is published for `**linux/amd64`** only (see `environment.yaml`, `docker/`).
Apple Silicon (**arm64**) hosts resolve the wrong manifest by default — always
pull with an explicit CPU architecture:

```bash
docker pull --platform linux/amd64 danielnguyener/riboflow:latest
```

Then start an interactive shell with your clone of this repo mounted to
`/work`:

```bash
cd /path/to/riboflow_genome
docker run --platform linux/amd64 --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/work -w /work \
  danielnguyener/riboflow:latest bash
```

Inside the container, the image exposes the `ribo_genome` conda env on `PATH`;
`nextflow` **19.04.1** is installed there (matching the pipeline requirement).
If you open a fresh interactive shell where that `PATH` is not applied, run
`conda activate ribo_genome` before invoking Nextflow.

```bash
nextflow run RiboFlow.groovy -params-file example.yaml -c configs/local.config
```

`**--platform` placement:** Use it on both `pull` and `run`, **before** the
image name. Anything after `docker run … IMAGE` is interpreted as the
command executed *inside* the container — e.g.
`docker run IMAGE --platform …` fails because `--platform` is not a Docker
flag there.

Mounts for data outside the repo: add more `-v /host/dir:/mnt/dir`, then point
paths in your params YAML at `/mnt/dir/…`.

## Test Run

See [Four ways to run](#four-ways-to-run-local-or-hpc-conda-or-container) for when to choose **conda vs Docker vs Apptainer**. A minimal first run:

```bash
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/DanielNguyener/riboflow_genome.git
cd riboflow_genome

# Sample FASTQs (ribo + RNA-seq).
git clone https://github.com/DanielNguyener/rf_sample_data_genome.git

# NOTE: rf_sample_data includes a chrM genome index.
# See [Building the STAR genome index](#building-the-star-genome-index)
# (--sjdbOverhang 28 for ribo-seq) or use example_local.yaml + the chrM helper
# under references_for_riboflow/genome/. Bowtie paths are analogous under
# input.reference.filter.
```

Then either activate the conda env on a Linux host:

```bash
conda activate ribo_genome
nextflow run RiboFlow.groovy -params-file example.yaml
```

…or enter the apptainer container (typical **HPC clusters**):

```bash
# on TACC: module load tacc-apptainer/1.4.1
apptainer pull docker://danielnguyener/riboflow:latest
apptainer shell riboflow_latest.sif
nextflow run RiboFlow.groovy -params-file example.yaml
```

…or use **Docker Desktop** on macOS, Windows, or Linux workstations (Apple
Silicon hosts need `--platform linux/amd64`; see [Docker Desktop](#docker-desktop-local-workstations)):

```bash
docker pull --platform linux/amd64 danielnguyener/riboflow:latest
docker run --platform linux/amd64 --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/work -w /work \
  danielnguyener/riboflow:latest bash
```

Inside that shell:

```bash
nextflow run RiboFlow.groovy -params-file example.yaml -c configs/local.config
```

`example.yaml` is configured to write its outputs to `test_output/` and its
intermediates to `test_intermediates/` (set via `output.output.base` and
`output.intermediates.base`). Real project params files typically set these
back to `output` and `intermediates`. See the
[Output section](#output) for details.

## Running on an HPC Cluster (Apptainer / Singularity)

For HPC environments where Docker isn't available (TACC, most academic
clusters), run the pipeline **from inside an Apptainer container shell**
rather than letting Nextflow launch a fresh `apptainer exec` per task — that’s
the **HPC + container** cell in [Four ways to run](#four-ways-to-run-local-or-hpc-conda-or-container).
The same Linux nodes can alternatively use your own `**ribo_genome`** conda
install (**HPC + conda** in that table): tools on `PATH`, then identical
`nextflow run …`. Apptainer remains preferable on Lustre when you hit
squashfuse instability (explained below).

### Why an interactive shell rather than `-profile singularity_`*?

Nextflow 19.04.1's `singularity` scope wraps every process in a separate
`singularity exec` invocation. On Lustre-backed shared filesystems we hit two
failure modes with long-running tasks:

1. `**Transport endpoint is not connected` (ENOTCONN)** from `/usr/bin/<tool>`
  — squashfuse loses its backing under concurrent I/O on Lustre.
2. **Silent PATH degradation** — bash mid-task reports
  `awk: command not found` (and friends) because squashfuse `readdir` on
   `/usr/bin` starts returning nothing partway through the run.

Launching one `apptainer shell` and running `nextflow run ...` *inside* it
gives the whole pipeline a single, stable squashfuse mount held by your
interactive shell, sidestepping both failure modes.

### One-time setup

```bash
apptainer pull docker://danielnguyener/riboflow:latest
```

This writes `riboflow_latest.sif` to the current directory. On TACC LS6
the system `apptainer.conf` auto-binds `/scratch`, `/work`, and `/home1`,
so no `--bind` flag is needed at pull or shell time. On other clusters
you may need to add `--bind /path/to/data:/path/to/data` to the
`apptainer shell` invocation below.

### Per-run workflow

```bash
apptainer shell riboflow_latest.sif

# Now inside the container — every tool (samtools, STAR, bamCoverage, ...)
# is on PATH from the container's conda env.
cd /path/to/your_run_dir
nextflow run /path/to/riboflow_genome/RiboFlow.groovy \
    -params-file /path/to/your_params.yaml
```

Inside the shell, Nextflow runs with the default profile (no `singularity`
scope enabled), so every process is just a regular subprocess of your shell.
That's the whole trick.

### Non-interactive batch (`sbatch` / `apptainer exec`)

If you need to launch the pipeline non-interactively (e.g. from an SBATCH
script), wrap the whole job in a single `apptainer exec` so the same shell
holds the squashfuse mount for the lifetime of the run:

```bash
#!/usr/bin/env bash
#SBATCH -J riboflow
#SBATCH -t 24:00:00
module load tacc-apptainer/1.4.1
apptainer exec --bind /scratch/ riboflow_latest.sif bash -c '
    cd /path/to/your_run_dir
    nextflow run /path/to/riboflow_genome/RiboFlow.groovy \
        -params-file /path/to/your_params.yaml
'
```

Avoid spawning a fresh `apptainer exec` per Nextflow task (which is what the
old `singularity_cluster` profile did) — that's the failure mode this
recipe is designed to sidestep.

## Output

The base output and intermediates directories are set in your params file:

```yaml
output:
   individual_lane_directory: 'individual'
   merged_lane_directory: 'merged'
   intermediates:
      base: 'test_intermediates'   # → $NF_RUN_DIR/test_intermediates/
   output:
      base: 'test_output'          # → $NF_RUN_DIR/test_output/
```

`example.yaml` writes to `test_output/` / `test_intermediates/`. Production
runs typically use `output` / `intermediates`. The trees below use `<out>` and
`<inter>` for whatever you configure.

### Output directory (`<out>/`)

The exact files depend on `dedup_method` and which optional features are
enabled. The two most common configurations are shown below.

#### `dedup_method: "umicollapse"` with `do_rnaseq: true`, `do_strand_split: true`

```
<out>/
├── alignments/
│   ├── ribo/
│   │   ├── individual/
│   │   │   ├── <sample>.<lane>.genome.qpass.bed
│   │   │   ├── <sample>.<lane>.genome.post_dedup.bed
│   │   │   ├── <sample>.<lane>.post_dedup.bam
│   │   │   └── <sample>.<lane>.post_dedup.bam.bai
│   │   ├── merged/
│   │   │   ├── <sample>.dedup.bam
│   │   │   ├── <sample>.dedup.bam.bai
│   │   │   ├── <sample>.genome.post_dedup.bed
│   │   │   ├── <sample>.genome.qpass.merged.bam
│   │   │   └── <sample>.genome.qpass.merged.bam.bai
│   │   └── stranded/                            # only if do_strand_split: true
│   │       ├── <sample>.ribo.plus.bam
│   │       ├── <sample>.ribo.plus.bam.bai
│   │       ├── <sample>.ribo.plus.bed
│   │       ├── <sample>.ribo.minus.bam
│   │       ├── <sample>.ribo.minus.bam.bai
│   │       └── <sample>.ribo.minus.bed
│   └── rnaseq/                                  # only if do_rnaseq: true
│       ├── individual/
│       │   └── <sample>.<lane>.rnaseq_genome.qpass.bed
│       └── merged/
│           ├── <sample>.rnaseq_genome.qpass.bed
│           ├── <sample>.rnaseq_genome.qpass.merged.bam
│           └── <sample>.rnaseq_genome.qpass.merged.bam.bai
├── bigwigs/
│   ├── ribo/
│   │   ├── <sample>.ribo.plus.bigWig
│   │   └── <sample>.ribo.minus.bigWig
│   └── rnaseq/                                  # only if do_rnaseq: true
│       └── <sample>.rnaseq.bigWig
├── rnaseq/                                      # only if do_rnaseq: true
│   └── stats/
│       ├── rnaseq_stats.csv
│       └── rnaseq_individual_stats.csv
└── stats/
    ├── stats.csv
    ├── individual_stats.csv
    └── index_fastq_correspondence.txt
```

#### `dedup_method: "position"` with `do_rnaseq: true`, `do_strand_split: true`

The ribo-seq individual directory contains BEDs only (no per-lane BAM; the
position deduplicator operates on a merged BED). The merged directory gains
both the post-dedup BAM and BED. The stranded and bigwig outputs are
identical in shape to the umicollapse case.

```
<out>/
├── alignments/
│   ├── ribo/
│   │   ├── individual/
│   │   │   ├── <sample>.<lane>.genome.qpass.bed
│   │   │   └── <sample>.<lane>.genome.post_dedup.bed
│   │   ├── merged/
│   │   │   ├── <sample>.genome.post_dedup.bed
│   │   │   ├── <sample>.genome.qpass.merged.bam
│   │   │   ├── <sample>.genome.qpass.merged.bam.bai
│   │   │   ├── <sample>.post_dedup.bam
│   │   │   └── <sample>.post_dedup.bam.bai
│   │   └── stranded/                            # only if do_strand_split: true
│   │       ├── <sample>.ribo.plus.bam
│   │       ├── <sample>.ribo.plus.bam.bai
│   │       ├── <sample>.ribo.plus.bed
│   │       ├── <sample>.ribo.minus.bam
│   │       ├── <sample>.ribo.minus.bam.bai
│   │       └── <sample>.ribo.minus.bed
│   └── rnaseq/                                  # only if do_rnaseq: true
│       ├── individual/
│       │   ├── <sample>.<lane>.rnaseq_genome.qpass.bed
│       │   └── <sample>.<lane>.rnaseq_genome.post_dedup.bed
│       └── merged/
│           ├── <sample>.rnaseq_genome.post_dedup.bam
│           ├── <sample>.rnaseq_genome.post_dedup.bam.bai
│           ├── <sample>.rnaseq_genome.post_dedup.bed
│           ├── <sample>.rnaseq_genome.qpass.merged.bam
│           └── <sample>.rnaseq_genome.qpass.merged.bam.bai
├── bigwigs/
│   ├── ribo/
│   │   ├── <sample>.ribo.plus.bigWig
│   │   └── <sample>.ribo.minus.bigWig
│   └── rnaseq/                                  # only if do_rnaseq: true
│       └── <sample>.rnaseq.bigWig
├── rnaseq/                                      # only if do_rnaseq: true
│   └── stats/
│       ├── rnaseq_stats.csv
│       └── rnaseq_individual_stats.csv
└── stats/
    ├── stats.csv
    ├── individual_stats.csv
    └── index_fastq_correspondence.txt
```

### Intermediates directory (`<inter>/`)

All intermediate files are safe to delete; Nextflow's `storeDir` will
regenerate them on re-run. The layout mirrors the two dedup cases.

#### `dedup_method: "umicollapse"`

```
<inter>/
├── alignment_ribo/
│   ├── individual/
│   │   ├── <sample>.<lane>.dedup.{total,primary,secondary}.count
│   │   ├── <sample>.<lane>.genome.post_dedup.bed
│   │   ├── <sample>.<lane>.post_dedup.bam
│   │   └── <sample>.<lane>.post_dedup.bam.bai
│   └── merged/
│       ├── <sample>.dedup.bam
│       ├── <sample>.dedup.bam.bai
│       └── <sample>.merged_dedup.{total,primary,secondary}.count
├── bam_to_bed/
│   └── individual/
│       └── <sample>.<lane>.genome.qpass.bed
├── clip/
│   ├── <sample>.<lane>.clipped.fastq.gz
│   └── <sample>.<lane>.clipped.log
├── filter/
│   ├── <sample>.<lane>.aligned.filter.fastq.gz
│   ├── <sample>.<lane>.filter.bam
│   ├── <sample>.<lane>.filter.log
│   └── <sample>.<lane>.unaligned.filter.fastq.gz
├── genome_alignment/
│   ├── individual/
│   │   ├── <sample>.<lane>.genome_alignment.aligned.fastq.gz
│   │   ├── <sample>.<lane>.genome_alignment.bam
│   │   ├── <sample>.<lane>.genome_alignment.log
│   │   ├── <sample>.<lane>.genome_alignment.secondary.count
│   │   └── <sample>.<lane>.genome_alignment.unaligned.fastq.gz
│   └── merged/
│       ├── <sample>.genome.qpass.merged.bam
│       └── <sample>.genome.qpass.merged.bam.bai
├── quality_filter/
│   ├── <sample>.<lane>.genome_alignment.qpass.bam
│   └── <sample>.<lane>.qpass.{total,primary,secondary}.count
├── umi_tools/
│   └── merged/
│       ├── <sample>.<lane>.umi_extracted.fastq.gz
│       └── <sample>.<lane>.umi_extracted.log
└── rnaseq/                                      # only if do_rnaseq: true
    ├── bam_to_bed/{individual,merged}/
    ├── clip/
    ├── filter/
    ├── genome_alignment/{individual,merged}/
    ├── quality_filter/
    └── stats/genome/{individual,merged}/
```

#### `dedup_method: "position"`

Adds position-dedup bookkeeping files. Also shows the
`star.output_transcriptome_bam: true` transcriptome intermediates.

```
<inter>/
├── alignment_ribo/
│   ├── individual/
│   │   ├── <sample>.<lane>.dedup.{total,primary,secondary}.count
│   │   └── <sample>.<lane>.genome.post_dedup.bed
│   └── merged/
│       ├── <sample>.merged_dedup.{total,primary,secondary}.count
│       ├── <sample>.post_dedup.bam
│       └── <sample>.post_dedup.bam.bai
├── bam_to_bed/
│   ├── individual/
│   │   ├── <sample>.<lane>.genome.qpass.bed
│   │   └── <sample>.<lane>.genome.with_sample_index.bed
│   └── merged/
│       └── <sample>.genome.merged.pre_dedup.bed
├── clip/
├── filter/
├── genome_alignment/
│   ├── individual/
│   │   ├── <sample>.<lane>.genome_alignment.{bam,log,secondary.count,...}
│   │   └── <sample>.<lane>.transcriptome_alignment.bam  # if output_transcriptome_bam: true
│   └── merged/
├── quality_filter/
├── transcriptome_alignment/              # if star.output_transcriptome_bam: true
│   ├── individual/
│   │   ├── <sample>.<lane>.transcriptome_alignment.qpass.bam
│   │   ├── <sample>.<lane>.transcriptome_alignment.qpass.bam.bai
│   │   ├── <sample>.<lane>.transcriptome.post_dedup.bed
│   │   └── <sample>.<lane>.transcriptome.with_sample_index.bed
│   └── merged/
│       ├── <sample>.transcriptome.merged.pre_dedup.bed
│       ├── <sample>.transcriptome.post_dedup.bam
│       ├── <sample>.transcriptome.post_dedup.bam.bai
│       ├── <sample>.transcriptome.post_dedup.bed
│       ├── <sample>.transcriptome.qpass.merged.bam
│       └── <sample>.transcriptome.qpass.merged.bam.bai
└── rnaseq/                              # only if do_rnaseq: true
    ├── bam_to_bed/{individual,merged}/
    ├── clip/
    ├── filter/
    ├── genome_alignment/{individual,merged}/
    ├── quality_filter/
    └── stats/genome/{individual,merged}/
```

Ribo-seq bigwigs cover read 5' ends on the genome (no P-site correction;
this is a pure genome-alignment pipeline). RNA-seq bigwigs are unstranded
coverage.

### Stats CSV schema

`stats.csv` and `individual_stats.csv` are wide-format (one column per
sample/lane). Row labels follow this pattern, where `<step>` is one of
`genome`, `qpass`, or `dedup`:

- `total_reads`, `clipped_reads`, `filtered_out`, `filter_kept`,
`genome_aligned_once`, `genome_aligned_many`, `genome_unaligned`
- `<step>_primary_alignments`, `<step>_secondary_alignments`,
`<step>_total_alignments`, `<step>_pct_primary`
(within-step primary share = primary / total)
- `clipped_reads_%`, `filter_kept_%`, `filtered_out_%`,
`genome_primary_alignments_%`, `qpass_primary_alignments_%`,
`dedup_primary_alignments_%`
(reads-based retention from the previous step)

The percentage rows are computed by `scripts/stats_percentage.py`
during the merged-stats step.

### Compute requirements

Right-sizing depends on the genome and dataset, but a useful baseline
for the human genome:

- **STAR genome alignment**: ~30 GB resident for the loaded index plus
working memory for two-pass alignment. Roughly 10–15 min per ribo-seq
lane on 16 cores.
- **Bowtie2 rRNA filter**: ~2 GB RAM. The shipped pipeline caps
alignment threads at 16 and `samtools sort` threads at 8 per lane to
keep heap predictable.
- **umicollapse dedup**: ~32 GB JVM heap on CCLE-scale BAMs.
Memory-bound.
- **deepTools bamCoverage**: ~8 forked workers per bigwig process; cap
enforced in `RiboFlow.groovy` to avoid `fork()` ENOMEM.
- **Disk**: budget ~3–5× the input FASTQ size for `<inter>/`.

`configs/local.config` is sized for a 128-core / 256 GB TACC LS6 compute
node, with `maxForks` caps on the heavy processes
(`rnaseq_genome_alignment`, `genome_deduplicate_umicollapse`,
`transcriptome_deduplicate_umicollapse`, `genome_create_strand_specific_bigwigs`,
`rnaseq_create_bigwig`, `filter`, `rnaseq_filter`) so the node memory
budget can't be oversubscribed. For smaller workstations or other
clusters, copy `configs/local.config` to a sibling file, adjust
`cpus` / `memory` / `executor.memory` to match your node, and pass it
via `-c your_custom.config` on the `nextflow run` command line.

There is intentionally **no** `singularity_cluster` profile — per-task
`apptainer exec` is unreliable on Lustre. Use the
[apptainer-shell workflow](#running-on-an-hpc-cluster-apptainer--singularity)
instead.

## Building the STAR genome index

`input.reference.genome` must be a **directory** produced by STAR’s
`genomeGenerate` mode. The directory must contain at least `SA`, `SAindex`,
`Genome`, and `chrNameLength.txt`. Build the index with the **same STAR
major/minor** as the pipeline environment (`environment.yaml` pins
**STAR ≥ 2.7.10**).

### Full-genome build (GENCODE GRCh38 + annotation)

Decompress the FASTA and GTF (STAR expects plain-text inputs for
`--genomeFastaFiles` and `--sjdbGTFfile`), then run:

```bash
STAR --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir /path/to/STAR_GRCh38_index \
  --genomeFastaFiles /path/to/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /path/to/gencode.v48.annotation.gtf \
  --sjdbOverhang 28
```

Point `input.reference.genome` at `/path/to/STAR_GRCh38_index` (the directory
itself, not a file inside it).

### `--sjdbOverhang` for ribosome profiling

When you pass `--sjdbGTFfile`, STAR embeds annotated splice junctions into the
genome index. `**--sjdbOverhang` must be set to the maximum read length you will
align, minus 1**, taken over **every** library that shares this index (STAR
fixes the value at build time—there is no runtime override).

Ribosome-profiling footprints are short after adapter / quality trimming—often
on the order of **~26–34 nt**. **This project standardises on
`--sjdbOverhang 28` for ribo-seq–centric indices**, i.e. tuned for a **29 nt**
maximum read length in the footprint library. If your trimmed reads routinely
reach **32 nt**, rebuild with `**--sjdbOverhang 31`**; if you run both
ribo-seq and longer RNA-seq through the same STAR index, use `**(longest read in either assay) − 1`**.

If you omit `--sjdbGTFfile`, STAR builds a smaller index without annotated
junctions (usually a poor fit for spliced mammalian nuclear genes).

### Tiny index for local smoke tests

`references_for_riboflow/genome/scripts/build_star_mini_index.sh` wraps the same
`genomeGenerate` call for a **mitochondrial-only** (`chrM`, ~16 kb) slice,
defaulting to `**SJDB_OVERHANG=28`**. See `references_for_riboflow/genome/README.md`
and `example_local.yaml`. Nuclear reads align mostly as unmapped against that
mini index—use it to validate wiring, not biology.

### Fork `rf_sample_data` with your STAR genome

Upstream `[rf_sample_data](https://github.com/ribosomeprofiling/rf_sample_data)`
checks `mock_hg38` payloads into `**genome/`**, which targets **HISAT-style**
mock workflows—not **STAR**. To ship the upstream FASTQs, filter indices, and the
like unchanged while layering the chrM STAR directory you generated above, fork
`[ribosomeprofiling/rf_sample_data](https://github.com/ribosomeprofiling/rf_sample_data)`
on GitHub, clone **your fork** locally, and run:

```bash
./scripts/rf_sample_data_fork_replace_star_genome.sh /path/to/your_rf_sample_data_clone
```

By default this copies `references_for_riboflow/genome/STAR_index_chrM/` into the
fork’s `**genome/**` (override with `**STAR_DIR**`), writes `genome/README.md`,
and stages `git add genome`. Inspect `git status`, commit, and push. Then point params
YAML at `**./rf_sample_data/genome**` when the fork sits beside `**riboflow_genome**`.
The helper also supports cloning upstream automatically via `**BOOTSTRAP_CLONE_PARENT=…**`
(use only if GitHub/network can handle the transcriptome-heavy checkout).

## RiboFlow on Your Data

1. **Gather inputs.** You need:
  - Gzipped ribo-seq FASTQs (single-end).
  - (Optional) Gzipped RNA-seq FASTQs (single-end strings or
  paired-end `[R1, R2]` lists; sample names must match the ribo-seq
  names).
  - A Bowtie2 rRNA / contaminant filter index. The upstream
  [references_for_riboflow](https://github.com/ribosomeprofiling/references_for_riboflow)
  repo ships working indices for several organisms.
  - A STAR genome index **directory** (see [Building the STAR genome
  index](#building-the-star-genome-index); built with the STAR version
  pinned in `environment.yaml`; the directory must contain `SA`,
  `SAindex`, `Genome`, `chrNameLength.txt`).
2. **Copy `example.yaml` to `project.yaml`** and edit:
  - `input.reference.filter`, `input.reference.genome` — paths to the
   two reference resources above.
  - `input.fastq.<sample>` — one list per sample, each entry being a
  single-end string or a `[R1, R2]` paired-end list.
  - `rnaseq.fastq.<sample>` — same shape (or set `do_rnaseq: false`).
  - `dedup_method` (ribo) and `rnaseq.dedup_method` (RNA-seq).
  - `output.output.base` / `output.intermediates.base` — typically set
  these to `output` and `intermediates` for production runs (the
  shipped `test_output` / `test_intermediates` are deliberately
  namespaced so the sample run doesn't collide with real projects).
  - Adapter and alignment arguments under `clip_arguments`, `star.*`,
  and `alignment_arguments`. The shipped defaults are sensible for
  TruSeq small-RNA libraries with the `AAAAAAAAAACAAAAAAAAAA`
  poly-A adapter.
3. **Pick an execution environment.** Only one Nextflow profile (the
  `standard` default, `configs/local.config`) is shipped. Bring tools onto
   `PATH` using any **one** of the four patterns in
   [Four ways to run](#four-ways-to-run-local-or-hpc-conda-or-container):
  - **Local Linux:** `conda activate ribo_genome`.
  - **Local macOS/Windows:** `docker run --platform linux/amd64 … bash`
  then run Nextflow inside ([Docker Desktop](#docker-desktop-local-workstations)).
  - **HPC Linux (e.g. TACC):** either `conda activate ribo_genome` on the node,
  `**or`** `apptainer shell riboflow_latest.sif` (often preferable on Lustre;
  [HPC section](#running-on-an-hpc-cluster-apptainer--singularity)).
4. **Run.**
  ```bash
   nextflow run RiboFlow.groovy -params-file project.yaml
  ```
   Re-runs with `-resume` re-use cached steps, so you can iterate on
   downstream params (dedup method, qpass cutoffs, etc.) without
   re-aligning.

## Working with Unique Molecular Identifiers

Unique Molecular Identifiers (UMIs) are short random barcodes ligated to
each molecule before PCR. Reads sharing alignment position **and** UMI
are PCR duplicates and can be collapsed for more accurate quantification.

riboflow_genome deduplicates ribo-seq UMIs using
[umicollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse), a fast
Java-based tool that operates directly on coordinate-sorted BAMs. UMI
*extraction* (peeling the UMI off the 5′ end of the read into the FASTQ
header) is still handled by
[umi_tools extract](https://umi-tools.readthedocs.io/en/latest/), so two
params control the UMI flow:


| Param                         | Purpose                                                                                                                       |
| ----------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `umi_tools_extract_arguments` | Passed verbatim to `umi_tools extract`. Defines where the UMI lives on the read and how much surrounding sequence to discard. |
| `umicollapse_arguments`       | Extra flags forwarded to `umicollapse` on top of the core flags set in `RiboFlow.groovy`. Usually left as `""`.               |


### Enabling UMI dedup

In your params file:

```yaml
dedup_method: "umicollapse"

# 12 nt UMI at the 5' end, followed by 4 nt of spacer to discard.
umi_tools_extract_arguments: "-p \"^(?P<umi_1>.{12})(?P<discard_1>.{4}).+$\" --extract-method=regex"
umicollapse_arguments: ""
```

The shipped `example.yaml` already uses this layout (12 nt UMI + 4 nt
spacer, AAAAAAAAAACAAAAAAAAAA 3′ adapter) and exercises the full
umicollapse path — no separate `project_umi.yaml` is needed. To run it:

```bash
conda activate ribo_genome
nextflow run RiboFlow.groovy -params-file example.yaml
```

### Alternative ribo-seq dedup methods

`dedup_method` accepts three values:

- `"umicollapse"` — UMI-aware dedup. Requires the two params above.
- `"position"` — coordinate-only dedup (RFC `dedup`). Use when reads
have no UMI; alignments at the same start/strand are collapsed.
- `"none"` — skip dedup entirely. Bigwigs and final BEDs are then built
from the qpass BAM instead of a post-dedup BAM.

### UMI support for RNA-seq

UMIs are **not** supported on the RNA-seq side. RNA-seq lanes either
skip dedup (`rnaseq.dedup_method: "none"`) or use coordinate-based
dedup (`rnaseq.dedup_method: "position"`).

## A Note on References

**riboflow_genome** is designed to work with **genomic references** (STAR index +
Bowtie2 rRNA filter). This version is **NOT** configured for transcriptome alignment.

For transcriptome alignment, please use the original
[RiboFlow](https://github.com/ribosomeprofiling/riboflow).

## Advanced Features

### Transcriptome BAM output (`star.output_transcriptome_bam`)

Setting `star.output_transcriptome_bam: true` in your params file instructs STAR
to emit a second BAM in **transcriptome coordinates** (via `--quantMode
TranscriptomeSAM`) alongside the usual genome-coordinate BAM.

When this flag is `true` and `dedup_method` is `"umicollapse"` or `"position"`,
the pipeline runs the transcriptome BAM through the same quality-pass filter and
deduplication steps as the genome BAM. The deduped transcriptome-coordinate BAMs
land under:

```
<inter>/transcriptome_alignment/individual/
```

This is useful when you want to pass deduplicated, transcriptome-projected reads
to a downstream quantification tool (e.g. Salmon in alignment mode). When
`dedup_method: "none"`, the transcriptome BAMs are generated by STAR but not
further processed.

Set `star.output_transcriptome_bam: false` (the default) to skip this step
entirely and save disk space.

### Strand-split BAMs (`do_strand_split`)

Setting `do_strand_split: true` splits the merged, post-dedup ribo-seq BAM into
separate plus- and minus-strand BAMs using `samtools view` flag masks. The
split files are published under:

```
<out>/alignments/ribo/stranded/
  <sample>.ribo.plus.bam   (plus-strand reads)
  <sample>.ribo.plus.bam.bai
  <sample>.ribo.minus.bam  (minus-strand reads)
  <sample>.ribo.minus.bam.bai
```

The flag logic mirrors deepTools' `--filterRNAstrand` convention and is
strandedness-aware: for forward-stranded libraries (the default for
ribosome-profiling), plus-strand reads have flag 16 **not** set (`-F 2064`) and
minus-strand reads have flag 16 set (`-f 16 -F 2048`). The strandedness is
inferred automatically from your library layout.

These strand-split BAMs are independent of the bigWigs produced by the pipeline —
they give you per-strand BAM files you can load directly in a genome browser or
pass to another tool. `do_strand_split` defaults to `false`.

### Pairing Ribo-seq with RNA-seq

Set `do_rnaseq: true` in your params file and provide RNA-seq FASTQs under the
`rnaseq.fastq.<sample>` keys. Sample names must match the ribo-seq sample names.
RNA-seq lanes can be either single-end (string) or paired-end (`[R1, R2]` list).
See `example.yaml` for the full schema.

### Tuning quality-pass filters

`mapping_quality_cutoff` (ribo) and `rnaseq.mapping_quality_cutoff` set
the MAPQ threshold for the qpass step. `ribo_filter_flags` and
`rnaseq.filter_flags` are SAM flag masks passed to `samtools view -F`.

`example.yaml` ships with the unique-only configuration on both sides
(`ribo_filter_flags: 2308`, `rnaseq.filter_flags: 2308` — that's
4 unmapped + 256 secondary + 2048 supplementary). To let ribo-seq
multimappers contribute to bigwigs, switch the ribo mask to `2052`
(drop only 4 + 2048) and lower `mapping_quality_cutoff` to `0`. The
shipped `example.yaml` leaves the alternative commented out so you can
toggle it inline.

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

## [Frequently Asked Questions](FAQ.md)

## [Release Notes](CHANGELOG.md)

