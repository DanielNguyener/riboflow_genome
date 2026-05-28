# Channel contracts

## Canonical meta map

- **Per-lane** (one sequencing lane): `meta = [ id: <sample>, lane: <int>, strand: 'F' ]`
- **Per-sample** (after lane merge): `meta = [ id: <sample>, strand: 'F' ]`

Ribo-seq lanes are single-end (one FASTQ per lane). `strand` is always `'F'`
(forward-stranded) in this pipeline; it is carried in meta so the bigwig and
stranded-split steps need no separate strand channel. `type`/RNA-seq variants are
reserved for the deferred RNA-seq stage.

Per-lane channels are keyed (as element 0) by the **whole meta map**, so DSL2
`join` matches lanes when `[id,lane,strand]` are equal. Per-sample regrouping uses
`meta.id`: `ch.map{ meta,x -> [meta.id, x] }.groupTuple().map{ id, xs -> [[id:id, strand:'F'], xs] }`.

## Subworkflow I/O

### PREPROCESS(ch_reads, ch_filter_index)
- in `ch_reads`: `[ meta(lane), fastq ]`
- in `ch_filter_index` (value): `[ index_base, [bowtie2_index_files] ]`
- out `reads_for_genome`: `[ meta(lane), unaligned_fastq ]`
- out `clip_log`, `filter_log`: `[ meta(lane), log ]`

### GENOME_ALIGN(ch_reads_for_genome, ch_genome_index)
- in `ch_genome_index` (value): STAR index dir
- out `transcriptome_bam`: `[ meta(lane), txbam ]` (only if `star.output_transcriptome_bam`)
- out `genome_log`, `secondary_count`: `[ meta(lane), file ]`
- out `qpass_counts`, `individual_dedup_counts`: `[ meta(lane), total, primary, secondary ]`
- out `merged_dedup_counts`: `[ meta(sample), total, primary, secondary ]` (empty when dedup none)

### STAR_TRANSCRIPTOME_DEDUP(ch_tx_bam)  — optional, BAM/BED only
- in `ch_tx_bam`: `[ meta(lane), transcriptome_bam ]`
- out `bam`: `[ meta(sample), dedup_bam, bai ]`

### ALIGNMENT_STATS(clip_log, filter_log, genome_log, secondary_count, qpass_counts, individual_dedup_counts, merged_dedup_counts)
- out `individual`, `merged`: published `individual_stats.csv` / `stats.csv`

## Shared-module ext flags (set in conf/modules.config)

| ext flag | module(s) | effect |
|---|---|---|
| `presort` | SAMTOOLS_QPASS | pre-sort before filter (transcriptome BAM is QNAME-sorted) |
| `append_index` | BAM_TO_BED | append `sample.lane` as BED column 7 |
| `emit_counts` | SEPARATE_BED, RFC_EXTRACT_DEDUP_READS, UMICOLLAPSE_DEDUP | also emit total/primary/secondary counts |
| `emit_bed_counts` | SPLIT_DEDUP_BAM | also emit per-lane BED + counts |
| `use_merged_counts` | STATS_SUM | apply `rfc update-dedup-counts` (dedup ≠ none) |
| `prefix` | most | output basename (keeps DSL1 filenames) |
