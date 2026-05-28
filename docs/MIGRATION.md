# DSL1 → DSL2 migration map

Single source of truth for the migration. Each DSL1 process in `RiboFlow.groovy`
maps to a DSL2 module (`modules/local/*.nf`) or is deferred. A status box per row
lets the migration be executed/checked module-by-module without holding the
2653-line monolith in context.

Status: `[x]` migrated & stub-verified · `[~]` partial · `[ ]` deferred.

## Genome path (Phase 1)

| DSL1 process (RiboFlow.groovy line) | DSL2 module | Status |
|---|---|---|
| `write_fastq_correspondence` (182) | `write_correspondence.nf` | [x] |
| `raw_fastqc` (233) | `fastqc.nf` (alias RAW_FASTQC) | [x] |
| `clip` (257) | `cutadapt_clip.nf` | [x] |
| `extract_umi_via_umi_tools` (280) | `umitools_extract.nf` | [x] |
| `clipped_fastqc` (312) | `fastqc.nf` (alias CLIPPED_FASTQC) | [x] |
| `filter` (344) | `bowtie2_filter.nf` | [x] |
| `genome_alignment` (401) | `star_align.nf` | [x] |
| `genome_aligned_fastqc` (494) | `fastqc.nf` (alias GENOME_ALIGNED_FASTQC) | [x] |
| `genome_unaligned_fastqc` (523) | `fastqc.nf` (alias GENOME_UNALIGNED_FASTQC) | [x] |
| `genome_quality_filter` (576) | `samtools_qpass.nf` | [x] |
| `merge_genome_qpass_bam` (634) | `samtools_merge.nf` | [x] |
| `individual_genome_bam_to_bed` (983) | `bam_to_bed.nf` | [x] |
| `add_sample_index_col_to_genome_bed` (1012) | `add_sample_index_col.nf` | [x] |
| `merge_genome_bed_for_position_dedup` (1048) | `concat_sort_bed.nf` (alias MERGE_PRE_DEDUP_BED) | [x] |
| `genome_deduplicate_position` (1069) | `rfc_dedup.nf` | [x] |
| `separate_genome_bed_post_dedup` (1120) | `separate_bed.nf` (emit_counts) | [x] |
| `genome_convert_dedup_bed_to_bam_position` (1162) | `rfc_extract_dedup_reads.nf` (emit_counts) | [x] |
| `genome_deduplicate_umicollapse` (1208) | `umicollapse_dedup.nf` (emit_counts) | [x] |
| `split_genome_dedup_bam_to_individual` (1264) | `split_dedup_bam.nf` (emit_bed_counts) | [x] |
| `concat_genome_post_dedup_bed_umi` (1319) | `concat_sort_bed.nf` (alias CONCAT_POST_DEDUP_BED) | [x] |
| `concat_genome_qpass_bed_for_publish_none` (1356) | `concat_sort_bed.nf` (alias CONCAT_QPASS_BED_NONE) | [x] |
| `genome_create_strand_specific_bigwigs` (1391) | `deeptools_bamcoverage.nf` | [x] |
| `genome_split_stranded_bam` (1432) | `split_stranded_bam.nf` | [x] |
| `individual_genome_alignment_stats` (1528) | `stats_individual.nf` | [x] |
| `combine_individual_genome_alignment_stats` (1627) | `stats_combine.nf` (alias COMBINE_INDIVIDUAL) | [x] |
| `sum_individual_genome_alignment_stats` (1677) | `stats_sum.nf` | [x] |
| `combine_merged_genome_alignment_stats` (1718) | `stats_combine.nf` (alias COMBINE_MERGED) | [x] |
| `publish_stats` (1758) | `stats_publish.nf` | [x] |

## Transcriptome-projected BAM dedup (optional, `do_tx_dedup`)

All reuse the shared modules with transcriptome `ext` settings; BAM/BED only.

| DSL1 process (line) | DSL2 module | Status |
|---|---|---|
| `transcriptome_sort_and_filter` (675) | `samtools_qpass.nf` (presort) | [x] |
| `transcriptome_qpass_bam_to_bed` (715) | `bam_to_bed.nf` (append_index) | [x] |
| `merge_transcriptome_bed_for_position_dedup` (745) | `concat_sort_bed.nf` | [x] |
| `transcriptome_deduplicate_position` (764) | `rfc_dedup.nf` | [x] |
| `separate_transcriptome_bed_post_dedup` (801) | `separate_bed.nf` | [x] |
| `merge_transcriptome_qpass_bam_for_position` (830) | `samtools_merge.nf` | [x] |
| `transcriptome_convert_dedup_bed_to_bam_position` (859) | `rfc_extract_dedup_reads.nf` | [x] |
| `merge_transcriptome_bam_for_umicollapse` (896) | `samtools_merge.nf` | [x] |
| `transcriptome_deduplicate_umicollapse` (916) | `umicollapse_dedup.nf` | [x] |
| `split_transcriptome_dedup_bam_to_individual` (952) | `split_dedup_bam.nf` | [x] |

## Bowtie2 transcriptome alignment (new path, `transcriptome.run`)

| Component | DSL2 module / subworkflow | Status |
|---|---|---|
| bowtie2 transcriptome alignment | `bowtie2_transcriptome.nf` | [x] |
| transcriptome qpass filter | `samtools_qpass.nf` (ext.mapq) | [x] |
| transcriptome dedup fork (position/umicollapse/none) | `transcriptome_align.nf` (shared modules) | [x] |
| `ribopy create` / `.ribo` generation | `ribopy_create.nf` | [x] |
| transcriptome stats (per-lane + merged) | `tx_stats_individual.nf`, `transcriptome_stats.nf` | [x] |
| `genome.run` gate (genome-only / tx-only / both) | `nextflow.config` + `workflows/riboflow.nf` | [x] |

## Deferred (NOT migrated)

| Area | DSL1 location | Status |
|---|---|---|
| RNA-seq path (`rnaseq_*`) | `if (do_rnaseq)` 1781-2643 | [ ] |

## Key behavioural notes

- **umicollapse:** DSL1 ran a hand-shipped jar via a `java11` wrapper
  (`java11 -Xms512m -Xmx32g -Xss256m umicollapse.main.Main bam ...`). DSL2 uses
  the bioconda `umicollapse` CLI; JVM flags go via `_JAVA_OPTIONS` because the
  wrapper mis-parses `-Xss`. CLI flags are unchanged → output parity expected.
- **dedup_method** is global (not per-sample), resolved by `Utils.resolve_dedup_method`
  honouring the legacy `deduplicate` boolean. The fork lives in `genome_align.nf`.
- **stats join key** is the per-lane `meta`; dedup counts for `none` reuse qpass counts.

## Verification

- `nextflow config -profile test` parses. ✓
- Stub-run passes for dedup `none|position|umicollapse`, each ±`do_tx_dedup`. ✓
- Consolidated conda env (`environment.yaml`) + Docker build, umicollapse smoke
  test, and DSL1 parity diff must be run where conda/docker are available
  (see README / plan Verification §0, §5) — NOT possible in a no-daemon sandbox.
