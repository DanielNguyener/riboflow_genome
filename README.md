[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3376949.svg)](https://doi.org/10.5281/zenodo.3376949)

![RiboFlow](/docs/figures/riboflow_logo.jpg "RiboFlow Logo")

# RiboFlow_genome #

RiboFlow_genome is a forked version of the original [RiboFlow](https://github.com/ribosomeprofiling/riboflow) pipeline, customized for genomic alignment.
RiboFlow_genome belongs to a [software ecosystem](https://ribosomeprofiling.github.io/) desgined to work with ribosome profiling data.


![Overview](/docs/figures/ecosystem_overview.jpg "Ribo Ecosystem Overview")

## Contents

* [Installation](#installation)
* [Test Run](#test-run)
* [Running on an HPC Cluster (Apptainer / Singularity)](#running-on-an-hpc-cluster-apptainer--singularity)
* [Output](#output)
* [RiboFlow on Your Data](#riboflow-on-your-data)
* [UMIs](#working-with-unique-molecular-identifiers)
* [A Note on References](#a-note-on-references)
* [Advanced Features](#advanced-features)
* [FAQ](FAQ.md)
* [Changelog](CHANGELOG.md)

## Installation

### Requirements

* [Nextflow](https://www.nextflow.io/) **19.04.1** for running the DSL1 pipeline.
  Modern Nextflow (≥22) only parses DSL1 in config-check mode; actual runs need
  the 19.04.1 binary. Nextflow 19.04.1 requires Java 8 or 11.
* One of:
  * [Conda](https://conda.io/en/latest/miniconda.html) (recommended for local
    Linux workstations)
  * [Apptainer / Singularity](https://apptainer.org/) (HPC clusters — see the
    [HPC section below](#running-on-an-hpc-cluster-apptainer--singularity);
    this is also the recommended path on macOS, since the conda env is
    Linux-only)

### Conda option (recommended on Linux)

```bash
git clone https://github.com/DanielNguyener/riboflow_genome.git
conda env create -f riboflow_genome/environment.yaml
conda activate ribo_genome
```

The environment is named `ribo_genome` and ships every binary the
pipeline shells out to (STAR, bowtie2, samtools, cutadapt, umicollapse,
deeptools, bedtools, java11, …). Activate it before invoking Nextflow.

The conda env is **Linux-only** — several pinned packages don't build
on macOS. On a Mac (or any other host without Linux conda), use the
Apptainer workflow described in the
[HPC section](#running-on-an-hpc-cluster-apptainer--singularity).

### Which execution contexts are verified?

`tests/run_tests.sh` (Tier 3) runs the pipeline exactly once, in one of
two ways depending on where it's launched from:

| Context | When it runs | What the script does |
|---|---|---|
| `apptainer` | Launched from inside an `apptainer shell` (the container's `ribo_genome` env is already active) | Just runs `nextflow run RiboFlow.groovy -params-file example.yaml` |
| `conda` | Launched on a bare Linux host (no apptainer container around) | Sources `conda.sh` and `conda activate ribo_genome` itself, then runs nextflow |

Either path writes to `test_output/` / `test_intermediates/` under
`$NF_RUN_DIR` (the repo root by default). The conda path searches the
usual install locations for `conda.sh` (`$CONDA_EXE`, `conda info
--base`, `~/miniforge3`, `~/miniconda3`, `~/anaconda3`, `$WORK/miniforge3`);
if none are found, Tier 3 SKIPs with a clear message instead of failing
mid-run.

**No other execution contexts are tested** — in particular, the SLURM
executor and per-task `apptainer exec` (the deprecated
`singularity_cluster` profile) are deliberately not part of the test
suite. If you build something custom on top of those, you're on your own.

## Test Run

A first-time sanity check using the bundled `example.yaml`:

```bash
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/DanielNguyener/riboflow_genome.git
cd riboflow_genome

# Sample FASTQs (ribo + RNA-seq).
git clone https://github.com/ribosomeprofiling/rf_sample_data.git

# NOTE: rf_sample_data does NOT ship the STAR genome index (too large).
# Build one with `STAR --runMode genomeGenerate` against your reference
# of choice and point input.reference.genome at the resulting directory
# in your params file. The bowtie2 rRNA filter index path is set the
# same way under input.reference.filter.
```

Then either activate the conda env on a Linux host:

```bash
conda activate ribo_genome
nextflow run RiboFlow.groovy -params-file example.yaml
```

…or enter the apptainer container (recommended on the cluster and on
macOS — see the [HPC section](#running-on-an-hpc-cluster-apptainer--singularity)):

```bash
apptainer pull docker://danielnguyener/riboflow:latest
apptainer shell riboflow_latest.sif
nextflow run RiboFlow.groovy -params-file example.yaml
```

`example.yaml` is configured to write its outputs to `test_output/` and its
intermediates to `test_intermediates/` (set via `output.output.base` and
`output.intermediates.base`). Real project params files typically set these
back to `output` and `intermediates`. See the
[Output section](#output) for details.

## Running on an HPC Cluster (Apptainer / Singularity)

For HPC environments where Docker isn't available (TACC, most academic
clusters), run the pipeline **from inside an Apptainer container shell**
rather than letting Nextflow launch a fresh `apptainer exec` per task.

### Why an interactive shell rather than `-profile singularity_*`?

Nextflow 19.04.1's `singularity` scope wraps every process in a separate
`singularity exec` invocation. On Lustre-backed shared filesystems we hit two
failure modes with long-running tasks:

1. **`Transport endpoint is not connected` (ENOTCONN)** from `/usr/bin/<tool>`
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
module load tacc-apptainer
apptainer exec riboflow_latest.sif bash -c '
    cd /path/to/your_run_dir
    nextflow run /path/to/riboflow_genome/RiboFlow.groovy \
        -params-file /path/to/your_params.yaml
'
```

Avoid spawning a fresh `apptainer exec` per Nextflow task (which is what the
old `singularity_cluster` profile did) — that's the failure mode this
recipe is designed to sidestep.

## Output

The published directory names are controlled by the params file. The
relevant block looks like this in `example.yaml`:

```yaml
output:
   individual_lane_directory: 'individual'
   merged_lane_directory: 'merged'
   intermediates:
      base: 'test_intermediates'   # → $NF_RUN_DIR/test_intermediates/
   output:
      base: 'test_output'          # → $NF_RUN_DIR/test_output/
```

`example.yaml` writes to `test_output/` / `test_intermediates/` so the
shipped sample run doesn't collide with any real project tree. Typical
user params files set these to `output` / `intermediates`. The table below
uses `<out>` and `<inter>` to stand in for whatever you've configured.

| Path | Contents |
|---|---|
| `<out>/stats/stats.csv` | Per-sample alignment summary (one row per sample, including reads-based retention % and within-step primary-alignment %) |
| `<out>/stats/individual_stats.csv` | Per-lane alignment summary (same schema, one row per lane) |
| `<out>/bigwigs/ribo/*.ribo.{plus,minus}.bigWig` | Strand-specific ribo-seq bigWigs (built from the merged post-dedup BAM) |
| `<out>/bigwigs/rnaseq/*.rnaseq.bigWig` | RNA-seq coverage bigWigs (unstranded) |
| `<out>/alignments/ribo/individual/*.{bam,bed}` | Per-lane post-dedup ribo-seq alignments |
| `<out>/alignments/ribo/merged/*.{bam,bed}` | Merged-sample post-dedup ribo-seq alignments |
| `<out>/alignments/rnaseq/{individual,merged}/` | RNA-seq qpass / post-dedup alignments (depending on `rnaseq.dedup_method`) |
| `<out>/fastqc/` | FastQC reports (only if `do_fastqc: true`) |
| `<inter>/` | Cached working files (raw STAR BAMs, qpass BAMs, pre-dedup BEDs). Safe to delete; will be regenerated on re-run. |

Bigwigs and post-dedup BED/BAM artifacts are **only** emitted for the
final step of each branch — i.e. the post-deduplication outputs, or
the qpass outputs in the `dedup_method: "none"` branch. Intermediate
qpass / pre-dedup files stay under `<inter>/`.

Ribo-seq bigwigs cover read 5' ends on the genome (no P-site
correction; this is a pure genome-alignment pipeline). RNA-seq bigwigs
are unstranded coverage.

### Stats CSV schema

`stats.csv` and `individual_stats.csv` are wide-format (one column per
sample/lane). Row labels follow this pattern, where `<step>` is one of
`genome`, `qpass`, or `dedup`:

* `total_reads`, `clipped_reads`, `filtered_out`, `filter_kept`,
  `genome_aligned_once`, `genome_aligned_many`, `genome_unaligned`
* `<step>_primary_alignments`, `<step>_secondary_alignments`,
  `<step>_total_alignments`, `<step>_pct_primary`
  (within-step primary share = primary / total)
* `clipped_reads_%`, `filter_kept_%`, `filtered_out_%`,
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


## RiboFlow on Your Data

1. **Gather inputs.** You need:
   * Gzipped ribo-seq FASTQs (single-end).
   * (Optional) Gzipped RNA-seq FASTQs (single-end strings or
     paired-end `[R1, R2]` lists; sample names must match the ribo-seq
     names).
   * A Bowtie2 rRNA / contaminant filter index. The upstream
     [references_for_riboflow](https://github.com/ribosomeprofiling/references_for_riboflow)
     repo ships working indices for several organisms.
   * A STAR genome index **directory** (built with the STAR version
     pinned in `environment.yaml`; the directory must contain `SA`,
     `SAindex`, `Genome`, `chrNameLength.txt`).

2. **Copy `example.yaml` to `project.yaml`** and edit:
   * `input.reference.filter`, `input.reference.genome` — paths to the
     two reference resources above.
   * `input.fastq.<sample>` — one list per sample, each entry being a
     single-end string or a `[R1, R2]` paired-end list.
   * `rnaseq.fastq.<sample>` — same shape (or set `do_rnaseq: false`).
   * `dedup_method` (ribo) and `rnaseq.dedup_method` (RNA-seq).
   * `output.output.base` / `output.intermediates.base` — typically set
     these to `output` and `intermediates` for production runs (the
     shipped `test_output` / `test_intermediates` are deliberately
     namespaced so the sample run doesn't collide with real projects).
   * Adapter and alignment arguments under `clip_arguments`, `star.*`,
     and `alignment_arguments`. The shipped defaults are sensible for
     TruSeq small-RNA libraries with the `AAAAAAAAAACAAAAAAAAAA`
     poly-A adapter.

3. **Pick an execution environment.** Only one profile (the `standard`
   default, using `configs/local.config`) is shipped. Pick how you bring
   the binaries onto PATH:
   * **Conda** on a local Linux workstation: activate `ribo_genome`.
   * **Apptainer/Singularity** on an HPC cluster (or macOS): enter the
     container shell first. See the
     [HPC section](#running-on-an-hpc-cluster-apptainer--singularity).

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

| Param | Purpose |
|---|---|
| `umi_tools_extract_arguments` | Passed verbatim to `umi_tools extract`. Defines where the UMI lives on the read and how much surrounding sequence to discard. |
| `umicollapse_arguments` | Extra flags forwarded to `umicollapse` on top of the core flags set in `RiboFlow.groovy`. Usually left as `""`. |

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

* `"umicollapse"` — UMI-aware dedup. Requires the two params above.
* `"position"` — coordinate-only dedup (RFC `dedup`). Use when reads
  have no UMI; alignments at the same start/strand are collapsed.
* `"none"` — skip dedup entirely. Bigwigs and final BEDs are then built
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
