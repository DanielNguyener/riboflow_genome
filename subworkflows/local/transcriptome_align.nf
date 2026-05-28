// Bowtie2 transcriptome alignment → qpass → dedup → ribopy create (.ribo).
// Parallel branch to GENOME_ALIGN; both receive PREPROCESS.out.reads_for_genome.
// Gated by params.transcriptome.run in workflows/riboflow.nf.
//
// Dedup fork mirrors GENOME_ALIGN / STAR_TRANSCRIPTOME_DEDUP: position, umicollapse,
// or none. Shared modules are disambiguated via .*:TRANSCRIPTOME_ALIGN:MODULE
// selectors in conf/modules.config.
//
// Count simplification: transcriptome path uses total counts only (no primary/secondary).
// tx_stats_individual.nf only needs total qpass and total dedup — no merged-dedup
// override needed because per-lane totals sum correctly for transcriptome.

include { BOWTIE2_TRANSCRIPTOME }                          from '../../modules/local/bowtie2_transcriptome.nf'
include { RIBOPY_CREATE }                                   from '../../modules/local/ribopy_create.nf'
include { SAMTOOLS_QPASS }                                  from '../../modules/local/samtools_qpass.nf'
include { SAMTOOLS_MERGE }                                  from '../../modules/local/samtools_merge.nf'
include { BAM_TO_BED as TX_QPASS_BAM_TO_BED }              from '../../modules/local/bam_to_bed.nf'
include { BAM_TO_BED as TX_MERGED_QPASS_BED }              from '../../modules/local/bam_to_bed.nf'
include { BAM_TO_BED as TX_MERGED_DEDUP_BED }              from '../../modules/local/bam_to_bed.nf'
include { CONCAT_SORT_BED as MERGE_TX_PRE_DEDUP_BED }      from '../../modules/local/concat_sort_bed.nf'
include { RFC_DEDUP }                                       from '../../modules/local/rfc_dedup.nf'
include { SEPARATE_BED }                                    from '../../modules/local/separate_bed.nf'
include { RFC_EXTRACT_DEDUP_READS }                         from '../../modules/local/rfc_extract_dedup_reads.nf'
include { UMICOLLAPSE_DEDUP }                               from '../../modules/local/umicollapse_dedup.nf'
include { SPLIT_DEDUP_BAM }                                 from '../../modules/local/split_dedup_bam.nf'

workflow TRANSCRIPTOME_ALIGN {
    take:
    ch_reads        // [ meta(id,lane,strand), fastq ]
    ch_tx_index     // value: [index_base, index_files]
    ch_regions_bed  // value: path — annotation BED (5'UTR/CDS/3'UTR)
    ch_lengths_tsv  // value: path — transcript lengths TSV

    main:
    def dedup = Utils.resolve_dedup_method(params)

    BOWTIE2_TRANSCRIPTOME(ch_reads, ch_tx_index)
    ch_tx_bam = BOWTIE2_TRANSCRIPTOME.out.bam   // [ meta(lane), bam ]

    SAMTOOLS_QPASS(ch_tx_bam)
    ch_qpass_bam       = SAMTOOLS_QPASS.out.bam         // [ meta(lane), bam, bai ]
    ch_qpass_total     = SAMTOOLS_QPASS.out.total_count  // [ meta(lane), total ] (no primary/secondary)

    ch_lane_meta = ch_tx_bam.map { meta, bam -> [meta.id, meta] }

    // Merge per-lane qpass BAMs → per-sample.
    ch_merge_in = ch_qpass_bam
        .map { meta, bam, bai -> [meta.id, bam] }
        .groupTuple()
        .map { id, bams -> [[id: id, strand: 'F'], bams] }
    SAMTOOLS_MERGE(ch_merge_in)
    ch_merged_qpass_bam = SAMTOOLS_MERGE.out.bam  // [ smeta, bam, bai ]

    ch_ribo_bed             = Channel.empty()
    ch_individual_dedup_cnt = Channel.empty()   // [ meta, total ] — total only for transcriptome

    if (dedup == 'position') {
        // Per-lane qpass BED (with sample.lane appended in col 7) → merge → dedup.
        TX_QPASS_BAM_TO_BED(ch_qpass_bam.map { meta, bam, bai -> [meta, bam] })
        MERGE_TX_PRE_DEDUP_BED(
            TX_QPASS_BAM_TO_BED.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
        RFC_DEDUP(MERGE_TX_PRE_DEDUP_BED.out.bed)
        ch_ribo_bed = RFC_DEDUP.out.bed  // 7-col; RIBOPY_CREATE strips to 6 via cut

        // Fan merged post-dedup BED back to lanes → per-lane BED + counts (total only).
        ch_sep_in = ch_lane_meta
            .combine(RFC_DEDUP.out.bed.map { smeta, bed -> [smeta.id, bed] }, by: 0)
            .map { id, meta, bed -> [meta, bed] }
        SEPARATE_BED(ch_sep_in)
        ch_individual_dedup_cnt = SEPARATE_BED.out.counts.map { meta, t, p, s -> [meta, t] }

        // Extract reads matching dedup BED from merged qpass BAM → final BAM.
        ch_extract_in = RFC_DEDUP.out.bed
            .map { smeta, bed -> [smeta.id, bed] }
            .join(ch_merged_qpass_bam.map { smeta, bam, bai -> [smeta.id, bam] })
            .map { id, bed, bam -> [[id: id, strand: 'F'], bed, bam] }
        RFC_EXTRACT_DEDUP_READS(ch_extract_in)
    }
    else if (dedup == 'umicollapse') {
        UMICOLLAPSE_DEDUP(ch_merged_qpass_bam)

        // Split merged dedup BAM back to per-lane BAMs + counts (total only).
        ch_split_in = ch_lane_meta
            .combine(UMICOLLAPSE_DEDUP.out.bam.map { smeta, bam, bai -> [smeta.id, bam, bai] }, by: 0)
            .map { id, meta, bam, bai -> [meta, bam, bai] }
        SPLIT_DEDUP_BAM(ch_split_in)
        ch_individual_dedup_cnt = SPLIT_DEDUP_BAM.out.counts.map { meta, t, p, s -> [meta, t] }

        // Convert merged dedup BAM → BED for ribopy.
        TX_MERGED_DEDUP_BED(UMICOLLAPSE_DEDUP.out.bam.map { smeta, bam, bai -> [smeta, bam] })
        ch_ribo_bed = TX_MERGED_DEDUP_BED.out.bed
    }
    else {
        // none: convert merged qpass BAM directly → BED for ribopy.
        TX_MERGED_QPASS_BED(ch_merged_qpass_bam.map { smeta, bam, bai -> [smeta, bam] })
        ch_ribo_bed = TX_MERGED_QPASS_BED.out.bed
        ch_individual_dedup_cnt = ch_qpass_total  // already [ meta, total ]
    }

    RIBOPY_CREATE(ch_ribo_bed, ch_regions_bed, ch_lengths_tsv)

    emit:
    bowtie2_log             = BOWTIE2_TRANSCRIPTOME.out.log  // [ meta(lane), log ]
    qpass_total_count       = ch_qpass_total                 // [ meta(lane), total ]
    individual_dedup_counts = ch_individual_dedup_cnt        // [ meta, total ] — total only
    ribo                    = RIBOPY_CREATE.out.ribo          // [ smeta, ribo ]
}
