// RNA-seq genome alignment + qpass + dedup fork (none|position|umicollapse) +
// merge + bigwig. Parallel to genome_align.nf (ribo-seq); shared modules are
// disambiguated via .*:RNASEQ_GENOME_ALIGN:MODULE selectors in conf/modules.config.
// No stranded split for RNA-seq (can be added later).

include { STAR_ALIGN_RNASEQ }                          from '../../modules/local/star_align_rnaseq.nf'
include { SAMTOOLS_QPASS }                             from '../../modules/local/samtools_qpass.nf'
include { BAM_TO_BED }                                 from '../../modules/local/bam_to_bed.nf'
include { SAMTOOLS_MERGE }                             from '../../modules/local/samtools_merge.nf'
include { ADD_SAMPLE_INDEX_COL }                       from '../../modules/local/add_sample_index_col.nf'
include { CONCAT_SORT_BED as RNASEQ_MERGE_PRE_DEDUP_BED }   from '../../modules/local/concat_sort_bed.nf'
include { CONCAT_SORT_BED as RNASEQ_CONCAT_QPASS_BED_NONE } from '../../modules/local/concat_sort_bed.nf'
include { CONCAT_SORT_BED as RNASEQ_CONCAT_POST_DEDUP_BED } from '../../modules/local/concat_sort_bed.nf'
include { RFC_DEDUP }                                  from '../../modules/local/rfc_dedup.nf'
include { SEPARATE_BED }                               from '../../modules/local/separate_bed.nf'
include { RFC_EXTRACT_DEDUP_READS }                    from '../../modules/local/rfc_extract_dedup_reads.nf'
include { UMICOLLAPSE_DEDUP }                          from '../../modules/local/umicollapse_dedup.nf'
include { SPLIT_DEDUP_BAM }                            from '../../modules/local/split_dedup_bam.nf'
include { DEEPTOOLS_BAMCOVERAGE }                      from '../../modules/local/deeptools_bamcoverage.nf'

workflow RNASEQ_GENOME_ALIGN {
    take:
    ch_reads_for_genome   // [ meta(id,lane,strand,single_end), reads ]
    ch_genome_index       // value: genome dir

    main:
    def dedup       = Utils.resolve_rnaseq_dedup_method(params)
    def unique_only = ((params.rnaseq?.genome?.mapping_quality_cutoff ?: params.rnaseq?.mapping_quality_cutoff ?: 4) as int) > 0
    def zero_file   = file("${projectDir}/assets/zero.count")

    STAR_ALIGN_RNASEQ(ch_reads_for_genome, ch_genome_index)

    SAMTOOLS_QPASS(STAR_ALIGN_RNASEQ.out.bam)
    ch_qpass_bam        = SAMTOOLS_QPASS.out.bam
    ch_qpass_total      = SAMTOOLS_QPASS.out.total_count
    ch_qpass_primary    = SAMTOOLS_QPASS.out.primary_count
    ch_qpass_secondary  = SAMTOOLS_QPASS.out.secondary_count
    ch_qpass_unique_cnt = SAMTOOLS_QPASS.out.unique_count

    BAM_TO_BED(ch_qpass_bam.map { meta, bam, bai -> [meta, bam] })

    ch_merge_in = ch_qpass_bam
        .map { meta, bam, bai -> [meta.id, bam] }
        .groupTuple()
        .map { id, bams -> [[id: id, strand: 'F'], bams] }
    SAMTOOLS_MERGE(ch_merge_in)
    ch_merged_qpass_bam = SAMTOOLS_MERGE.out.bam

    ch_lane_meta = STAR_ALIGN_RNASEQ.out.bam.map { meta, bam -> [meta.id, meta] }

    ch_final_bam            = Channel.empty()
    ch_individual_dedup_cnt = Channel.empty()
    ch_merged_dedup_cnt     = Channel.empty()

    if (dedup == 'none') {
        ch_final_bam            = ch_merged_qpass_bam
        ch_individual_dedup_cnt = ch_qpass_total.join(ch_qpass_primary).join(ch_qpass_secondary).join(ch_qpass_unique_cnt)
        RNASEQ_CONCAT_QPASS_BED_NONE(
            BAM_TO_BED.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
    }
    else if (dedup == 'position') {
        ADD_SAMPLE_INDEX_COL(BAM_TO_BED.out.bed)
        RNASEQ_MERGE_PRE_DEDUP_BED(
            ADD_SAMPLE_INDEX_COL.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
        RFC_DEDUP(RNASEQ_MERGE_PRE_DEDUP_BED.out.bed)

        ch_sep_in = ch_lane_meta
            .combine(RFC_DEDUP.out.bed.map { smeta, bed -> [smeta.id, bed] }, by: 0)
            .map { id, meta, bed -> [meta, bed] }
        SEPARATE_BED(ch_sep_in)
        ch_individual_dedup_cnt = unique_only
            ? SEPARATE_BED.out.total_count.map { meta, t -> [meta, t, t, zero_file, t] }
            : SEPARATE_BED.out.total_count.join(SEPARATE_BED.out.detail_counts)

        ch_extract_in = RFC_DEDUP.out.bed
            .map { smeta, bed -> [smeta.id, bed] }
            .join(ch_merged_qpass_bam.map { smeta, bam, bai -> [smeta.id, bam] })
            .map { id, bed, bam -> [[id: id, strand: 'F'], bed, bam] }
        RFC_EXTRACT_DEDUP_READS(ch_extract_in)
        ch_final_bam        = RFC_EXTRACT_DEDUP_READS.out.bam
        ch_merged_dedup_cnt = unique_only
            ? RFC_EXTRACT_DEDUP_READS.out.total_count.map { meta, t -> [meta, t, t, zero_file, t] }
            : RFC_EXTRACT_DEDUP_READS.out.total_count.join(RFC_EXTRACT_DEDUP_READS.out.detail_counts)
    }
    else if (dedup == 'umicollapse') {
        UMICOLLAPSE_DEDUP(ch_merged_qpass_bam)
        ch_final_bam        = UMICOLLAPSE_DEDUP.out.bam
        ch_merged_dedup_cnt = unique_only
            ? UMICOLLAPSE_DEDUP.out.total_count.map { meta, t -> [meta, t, t, zero_file, t] }
            : UMICOLLAPSE_DEDUP.out.total_count.join(UMICOLLAPSE_DEDUP.out.detail_counts)

        ch_split_in = ch_lane_meta
            .combine(UMICOLLAPSE_DEDUP.out.bam.map { smeta, bam, bai -> [smeta.id, bam, bai] }, by: 0)
            .map { id, meta, bam, bai -> [meta, bam, bai] }
        SPLIT_DEDUP_BAM(ch_split_in)
        ch_individual_dedup_cnt = unique_only
            ? SPLIT_DEDUP_BAM.out.total_count.map { meta, t -> [meta, t, t, zero_file, t] }
            : SPLIT_DEDUP_BAM.out.total_count.join(SPLIT_DEDUP_BAM.out.detail_counts)

        RNASEQ_CONCAT_POST_DEDUP_BED(
            SPLIT_DEDUP_BAM.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
    }

    DEEPTOOLS_BAMCOVERAGE(ch_final_bam)

    emit:
    genome_log               = STAR_ALIGN_RNASEQ.out.log
    secondary_count          = STAR_ALIGN_RNASEQ.out.secondary_count
    qpass_total_count        = ch_qpass_total
    qpass_primary_count      = ch_qpass_primary
    qpass_secondary_count    = ch_qpass_secondary
    qpass_unique_count       = ch_qpass_unique_cnt
    individual_dedup_counts  = ch_individual_dedup_cnt
    merged_dedup_counts      = ch_merged_dedup_cnt
}
