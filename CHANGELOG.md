# Changelog

All notable changes to **riboflow_genome** are documented here. This pipeline is the
Nextflow **DSL2** rewrite of the original DSL1 `RiboFlow.groovy`; the `0.0.x` series
tracks that rewrite.

## 0.0.2

### Added
- **Bowtie2 transcriptome → `.ribo` path** (`transcriptome.run: true`): bowtie2
  alignment → qpass → dedup → `ribopy create` per sample, then `ribopy merge` into
  `all.ribo`.
- **RNA-seq path** (`do_rnaseq: true`): parallel clip → bowtie2 filter → STAR genome
  (ENCODE defaults) → dedup, plus an RNA-seq transcriptome path that embeds the RNA-seq
  BED into the matching `.ribo` via `ribopy rnaseq set`. Single-end and paired-end
  supported (PE + `umicollapse` is rejected up front).
- **Path gating**: `genome.run`, `transcriptome.run`, and `do_rnaseq` can be combined
  freely (genome-only, transcriptome-only, or both; with or without RNA-seq).
- **Two genome-stats modes**: unique-only (`genome.mapping_quality_cutoff: 255`) and
  multi-mapper (`0`), with per-stage unique/multi/secondary breakdowns in the latter.
- Single consolidated conda environment (`environment.yaml`): Nextflow ≥24, Java 17,
  `umicollapse` from bioconda. `docker/Dockerfile` builds purely from it.
- Profiles `local`, `conda`, `apptainer`, `docker`, `test`.
- Stub fixtures + GitHub Actions CI exercising every major path.

### Changed
- Migrated from the DSL1 monolith `RiboFlow.groovy` to an nf-core-style DSL2 layout
  (`main.nf` → `workflows/` → `subworkflows/local/` → `modules/local/`).
- `umicollapse` now comes from bioconda — the hand-shipped `umicollapse.jar` and
  `java11` wrapper are gone.

## 0.0.1

### Added
- STAR genome alignment path (ribo-seq): clip → bowtie2 rRNA filter → STAR → qpass →
  dedup (`umicollapse` / `position` / `none`) → merge → strand-specific bigWigs.
- UMI extraction (`umi_tools extract`) and UMI-aware deduplication (`umicollapse`).
- Optional STAR transcriptome-*projected* BAM dedup
  (`star.output_transcriptome_bam: true`, BAM/BED only — not `.ribo`).
- Per-lane and merged alignment statistics.

## 0.0.0

### Added
- Initial DSL2 scaffold forked from RiboFlow, targeting genome alignment.
