# FAQ

This is **riboflow_genome**, a fork of [RiboFlow][upstream] specialized for STAR-based **genome alignment** of ribosome profiling and RNA-seq data. If you want the original transcriptome-based pipeline (which produces `.ribo` files for downstream RiboPy/RiboR analysis), use the upstream RiboFlow.

[upstream]: https://github.com/ribosomeprofiling/riboflow

---

### 1. Does this pipeline produce `.ribo` files?

No. `.ribo` is a transcriptome-coordinate format and cannot represent genome-aligned reads. This pipeline emits **bigWigs**, **deduplicated BAM/BED files**, and **alignment statistics**. If you need a `.ribo` file for downstream RiboPy / RiboR analysis, use upstream RiboFlow.

### 2. What outputs does it produce?

```
output/
  alignments/
    ribo/{individual,merged}/      # dedup BAM + BED, per-lane and merged
    rnaseq/{individual,merged}/
  bigwigs/
    ribo/                          # *.ribo.plus.bigWig, *.ribo.minus.bigWig
    rnaseq/                        # *.rnaseq.bigWig
  stats/
    stats.csv                      # merged-sample alignment summary
    individual_stats.csv           # per-lane alignment summary
  rnaseq/stats/                    # RNA-seq stats CSVs
  fastqc/                          # raw + clipped FastQC reports (if do_fastqc: true)
```

Intermediate work files (raw STAR BAMs, qpass BAMs, pre-dedup BEDs) live under `intermediates/` and are not user-facing — but they are cached, so re-runs are fast.

### 3. Does it support UMIs?

Yes — set `dedup_method: "umicollapse"` in your config and provide `umi_tools_extract_arguments` describing your UMI layout. The pipeline runs `umi_tools extract` for clipping and `umicollapse` for deduplication. See `example.yaml` for a working invocation.

### 4. Does it support paired-end RNA-seq?

Yes. Each entry in `rnaseq.fastq.<sample>` can be either a single-end string or a `[R1, R2]` two-element list. Ribo-seq lanes are always single-end (RPFs are short fragments).

### 5. Can I run it without Docker?

Yes, with caveats. The conda environment in `environment.yaml` covers most dependencies. However:
- **`umicollapse.jar`** is not in bioconda. The Docker image ships it. For a conda-only install, download `umicollapse.jar` manually and ensure `java11` is on `PATH` so the dedup step can call it.
- The Docker image is the supported configuration; treat conda-only as best-effort.

### 6. What strandedness assumptions does the pipeline make?

- **Ribo-seq** is treated as forward-stranded (the sequenced read is sense to the RPF). Strand-specific bigwigs are emitted with `.plus.bigWig` and `.minus.bigWig` suffixes.
- **RNA-seq** is treated as unstranded — bigwigs are coverage-only. RNA-seq strand auto-detection (RSeQC `infer_experiment.py`) was removed in favor of treating bigwigs as a coverage diagnostic.

### 7. How are alignment statistics generated?

Per-lane stats are collected from:
- cutadapt log (input reads, clipped reads)
- bowtie2 filter log (rRNA filtered out, kept)
- STAR `Log.final.out`, parsed by `scripts/parse_star_log.py` (uniquely mapped, multi-loci, unmapped totals, mismatch rate)
- `samtools view -c` on the qpass BAM, cached as `<bam>.qpass.count`
- dedup count from `wc -l` on the post-dedup BED

These per-lane CSVs are aggregated by `rfc merge overall-stats` (from `RFCommands_genome`) into `output/stats/stats.csv` and `output/stats/individual_stats.csv`.

### 8. Where do I find references (rRNA filter, STAR genome index)?

You provide them. Set `input.reference.filter` to a bowtie2 index prefix for your contaminant filter (rRNA, mtRNA, ...) and `input.reference.genome` to a STAR index directory.

For human references, the upstream [references_for_riboflow](https://github.com/ribosomeprofiling/references_for_riboflow) repository has filter indices that work here. STAR genome indices need to be built once with the same STAR version pinned in `environment.yaml`.
