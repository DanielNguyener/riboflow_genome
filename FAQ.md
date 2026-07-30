# FAQ

This is **riboflow_genome** (Nextflow DSL2), a fork of [RiboFlow][upstream] that adds
STAR-based **genome alignment** of ribosome-profiling and RNA-seq data alongside the
classic bowtie2 **transcriptome → `.ribo`** path. The entry point is `main.nf`.

[upstream]: https://github.com/ribosomeprofiling/riboflow

---

### 1. Does this pipeline produce `.ribo` files?

Yes, when `transcriptome.run: true`. The transcriptome path aligns reads with bowtie2
against a transcriptome index and runs `ribopy create` to produce a per-sample `.ribo`,
then `ribopy merge` to produce `all.ribo`. The **genome** path (`genome.run: true`) is
separate and emits deduplicated BAM/BED, alignment stats, and — when `do_bigwig: true` —
bigWigs, but never `.ribo`. You can run either path or both.

### 2. What outputs does it produce?

Everything RNA-seq lives under a single `rnaseq/` subtree, mirroring the ribo-seq layout:

```
<out>/
  alignments/ribo/{individual,merged,stranded}/    # dedup BAM + BED
  bigwigs/ribo/                                    # *.ribo.plus/minus.bigWig — if do_bigwig: true
  fastqc/{raw,clipped,genome_*,transcriptome_*}/   # if do_fastqc: true
  ribo/                                            # *.ribo + all.ribo (if transcriptome.run)
  stats/genome/{stats.csv, individual_stats.csv}
  stats/transcriptome/{stats.csv, individual_stats.csv}   # if transcriptome.run: true
  rnaseq/                                          # if do_rnaseq: true
    alignments/genome/{individual,merged}/
    bigwigs/genome/                                # *.rna.plus/minus.bigWig — if do_bigwig: true
    fastqc/{raw,clipped}/                          # if do_fastqc: true
    stats/{genome,transcriptome}/{stats.csv, individual_stats.csv}
```

`do_bigwig` and `do_fastqc` both default to **false**.

Intermediate work files (raw STAR BAMs, qpass BAMs, pre-dedup BEDs) live under
`intermediates/` and are cached via `storeDir`, so re-runs are fast.

### 3. Does it support UMIs?

Yes — set `dedup_method: "umicollapse"` and provide `umi_tools_extract_arguments`
describing your UMI layout. The pipeline runs `umi_tools extract` to peel the UMI into
the read header and `umicollapse` to deduplicate. See `example_umi_uniq.yaml` for a
working invocation.

### 4. Does it support paired-end RNA-seq?

Yes, on the **genome** path. Each entry in `rnaseq.fastq.<sample>` can be a single-end
string or a `[R1, R2]` two-element list. Ribo-seq lanes are always single-end (RPFs are
short fragments). Two PE combinations are rejected at startup
(`workflows/riboflow.nf:158-171`):

- `rnaseq.dedup_method: "position"` — needs fragment/BEDPE dedup. Use `umicollapse`
  (UMIs are extracted from both mates and deduped with `umicollapse --paired`; counts are
  fragment-level) or `none`.
- The RNA-seq **transcriptome** path under PE — set `transcriptome.run: false`, or use
  single-end.

`example_rnaseq_pe.yaml` is a working PE example.

### 5. Can I run it without Docker?

Yes. `environment.yaml` is a single consolidated conda env (`ribo_genome`) that ships
every tool, including `umicollapse` from **bioconda** — there is no longer a
hand-shipped `umicollapse.jar` or `java11` wrapper. Use `-profile conda` (Nextflow
builds the env) or `-profile default` with the env already activated. The conda env is
Linux-only (several pinned packages don’t build on macOS); on macOS/Windows use the
Docker/Apptainer image.

### 6. How are alignment statistics generated?

Per-lane stats are collected from the cutadapt log (input/clipped reads), the bowtie2
filter log (rRNA filtered out/kept), STAR’s `Log.final.out` (parsed by
`rfc parse-star-log`), and the qpass/dedup count files. The stats subworkflow combines
the per-lane CSVs and computes percentage rows (via the `rfc` helper from
`RFCommands_genome`) into `stats/genome/stats.csv` and `individual_stats.csv`. There is
no standalone `scripts/` directory.

### 7. Where do I find references (rRNA filter, STAR genome index, transcriptome)?

You provide them. `input.reference.filter` is a bowtie2 contaminant-filter index prefix;
`input.reference.genome` is a STAR index directory (see the README for how to build
one). For the `.ribo` path, also provide `transcriptome`, `regions`, and
`transcript_lengths`. The upstream
[references_for_riboflow](https://github.com/ribosomeprofiling/references_for_riboflow)
repository has filter and transcriptome references that work here. For the STAR genome
index you have two options: point `input.reference.genome` at one you built yourself with
the STAR version pinned in `environment.yaml` (Mode A), or set
`input.reference.genome_fasta` + `gtf` and let the pipeline run `genomeGenerate` for you
(Mode B) — see "Building the STAR genome index" in the README.
