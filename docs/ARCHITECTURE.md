# Architecture (Phase 1 — genome path)

`main.nf` → `workflows/riboflow.nf` builds the per-lane input channel, runs file
existence checks, then calls the subworkflows below. Environment selection
(conda/docker/ambient) is profile-driven; per-module output locations and CLI
args come from `conf/modules.config`.

```mermaid
flowchart TD
    IN["fastq: sample → [lane files]\n(riboflow.nf input build)"] --> CLIP
    subgraph PRE["PREPROCESS"]
      CLIP["CUTADAPT_CLIP"] -->|umicollapse| UMI["UMITOOLS_EXTRACT"]
      CLIP -->|none/position| FILT
      UMI --> FILT["BOWTIE2_FILTER (rRNA/tRNA)\n→ unaligned reads"]
    end
    FILT --> STAR["STAR_ALIGN\ngenome BAM (+optional tx-projected BAM)"]
    subgraph GEN["GENOME_ALIGN"]
      STAR --> QP["SAMTOOLS_QPASS\nmapq + flag filter, counts"]
      QP --> B2B["BAM_TO_BED (per lane)"]
      QP --> MERGE["SAMTOOLS_MERGE (per sample)"]
      MERGE --> DEDUP{"dedup_method"}
      DEDUP -->|none| FINAL["final sample BAM"]
      DEDUP -->|position| POS["ADD_SAMPLE_INDEX_COL → MERGE_PRE_DEDUP_BED\n→ RFC_DEDUP → SEPARATE_BED\n→ RFC_EXTRACT_DEDUP_READS"]
      DEDUP -->|umicollapse| UC["UMICOLLAPSE_DEDUP → SPLIT_DEDUP_BAM"]
      POS --> FINAL
      UC --> FINAL
      FINAL --> BW["DEEPTOOLS_BAMCOVERAGE\nribo.plus/minus.bigWig"]
      FINAL --> SS["SPLIT_STRANDED_BAM (do_strand_split)"]
    end
    STAR -. "tx-projected BAM (do_tx_dedup)" .-> STX["STAR_TRANSCRIPTOME_DEDUP\nBAM/BED only — NO .ribo"]
    subgraph ST["ALIGNMENT_STATS"]
      QP --> SI["STATS_INDIVIDUAL (per lane)\nrfc parse-star-log + python"]
      SI --> CI["COMBINE_INDIVIDUAL"]
      SI --> SUM["STATS_SUM (per sample)"]
      SUM --> CM["COMBINE_MERGED"]
      CI --> PUB["STATS_PUBLISH\nstats.csv, individual_stats.csv"]
      CM --> PUB
    end
```

## Output layout (parity with DSL1)
- `output/alignments/ribo/{individual,merged}/` — post-dedup BAM/BED (or qpass when `none`)
- `output/alignments/ribo/stranded/` — stranded BAM/BED (when `do_strand_split`)
- `output/bigwigs/ribo/*.ribo.{plus,minus}.bigWig` (only when `do_bigwig`; default off)
- `output/stats/{stats.csv,individual_stats.csv,index_fastq_correspondence.txt}`
- `output/fastqc/...` (when `do_fastqc`)
- intermediates under `intermediates/` (incl. `transcriptome_alignment/` for the tx path)

## Storage notes
- Final BAM/BED (and the `.ribo` from `RIBOPY_RNASEQ_SET`) are cached in `intermediates/`
  via `storeDir` and exposed under `output/` via `publishDir mode: 'link'` — a **hard link**,
  so there is only **one physical copy** on disk (no duplication). This requires `output/` and
  `intermediates/` to live on the **same filesystem** (the default — they're siblings in the
  run dir). If you point their `base` paths at different mounts, switch those `publishDir`
  entries to `mode: 'symlink'` in `conf/modules.config`.
- bigWig generation (`DEEPTOOLS_BAMCOVERAGE`, genome ribo + rnaseq) is gated on `do_bigwig`
  (default `false`); set `do_bigwig: true` to produce them.

See `` for the process→module map and ``
for the meta map and channel shapes.
