// Genome alignment + qpass + dedup fork (none|position|umicollapse) + merge +
// bigwig + stranded split. (RiboFlow.groovy:401-1464)

include { STAR_ALIGN }                          from '../../modules/local/star_align.nf'
include { FASTQC as GENOME_ALIGNED_FASTQC }     from '../../modules/local/fastqc.nf'
include { FASTQC as GENOME_UNALIGNED_FASTQC }   from '../../modules/local/fastqc.nf'
include { SAMTOOLS_QPASS }                      from '../../modules/local/samtools_qpass.nf'
include { BAM_TO_BED }                          from '../../modules/local/bam_to_bed.nf'
include { SAMTOOLS_MERGE }                      from '../../modules/local/samtools_merge.nf'
include { ADD_SAMPLE_INDEX_COL }                from '../../modules/local/add_sample_index_col.nf'
include { CONCAT_SORT_BED as MERGE_PRE_DEDUP_BED }   from '../../modules/local/concat_sort_bed.nf'
include { CONCAT_SORT_BED as CONCAT_QPASS_BED_NONE } from '../../modules/local/concat_sort_bed.nf'
include { CONCAT_SORT_BED as CONCAT_POST_DEDUP_BED } from '../../modules/local/concat_sort_bed.nf'
include { RFC_DEDUP }                           from '../../modules/local/rfc_dedup.nf'
include { SEPARATE_BED }                        from '../../modules/local/separate_bed.nf'
include { RFC_EXTRACT_DEDUP_READS }             from '../../modules/local/rfc_extract_dedup_reads.nf'
include { UMICOLLAPSE_DEDUP }                   from '../../modules/local/umicollapse_dedup.nf'
include { SPLIT_DEDUP_BAM }                     from '../../modules/local/split_dedup_bam.nf'
include { DEEPTOOLS_BAMCOVERAGE }               from '../../modules/local/deeptools_bamcoverage.nf'
include { SPLIT_STRANDED_BAM }                  from '../../modules/local/split_stranded_bam.nf'

workflow GENOME_ALIGN {
    take:
    ch_reads_for_genome   // [ meta(id,lane,strand), fastq ]
    ch_genome_index       // value: genome dir

    main:
    def dedup = Utils.resolve_dedup_method(params)

    STAR_ALIGN(ch_reads_for_genome, ch_genome_index)
    GENOME_ALIGNED_FASTQC(STAR_ALIGN.out.aligned)
    GENOME_UNALIGNED_FASTQC(STAR_ALIGN.out.unaligned)

    SAMTOOLS_QPASS(STAR_ALIGN.out.bam)
    ch_qpass_bam    = SAMTOOLS_QPASS.out.bam            // [ meta, bam, bai ]
    ch_qpass_counts = SAMTOOLS_QPASS.out.counts         // [ meta, total, primary, secondary ]

    // Per-lane qpass BED (published as the final per-lane artifact when none).
    BAM_TO_BED(ch_qpass_bam.map { meta, bam, bai -> [meta, bam] })

    // Merge per-lane qpass BAMs → per-sample qpass BAM.
    ch_merge_in = ch_qpass_bam
        .map { meta, bam, bai -> [meta.id, bam] }
        .groupTuple()
        .map { id, bams -> [[id: id, strand: 'F'], bams] }
    SAMTOOLS_MERGE(ch_merge_in)
    ch_merged_qpass_bam = SAMTOOLS_MERGE.out.bam        // [ smeta, bam, bai ]

    // Per-lane meta keyed by sample id, for fanning sample-level results back out.
    ch_lane_meta = STAR_ALIGN.out.bam.map { meta, bam -> [meta.id, meta] }

    // Defaults; each dedup branch overrides.
    ch_final_bam            = Channel.empty()   // [ smeta, bam, bai ] → bigwig/strand
    ch_individual_dedup_cnt = Channel.empty()   // [ meta, total, primary, secondary ] per lane
    ch_merged_dedup_cnt     = Channel.empty()   // [ smeta, total, primary, secondary ] per sample

    if (dedup == 'none') {
        ch_final_bam            = ch_merged_qpass_bam
        ch_individual_dedup_cnt = ch_qpass_counts
        // Publish merged qpass BED (concat of per-lane qpass BEDs).
        CONCAT_QPASS_BED_NONE(
            BAM_TO_BED.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
    }
    else if (dedup == 'position') {
        // Tag per-lane BEDs, merge, rfc-dedup.
        ADD_SAMPLE_INDEX_COL(BAM_TO_BED.out.bed)
        MERGE_PRE_DEDUP_BED(
            ADD_SAMPLE_INDEX_COL.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
        RFC_DEDUP(MERGE_PRE_DEDUP_BED.out.bed)

        // Fan merged post-dedup BED back to each lane → per-lane BED + counts.
        ch_sep_in = ch_lane_meta
            .combine(RFC_DEDUP.out.bed.map { smeta, bed -> [smeta.id, bed] }, by: 0)
            .map { id, meta, bed -> [meta, bed] }
        SEPARATE_BED(ch_sep_in)
        ch_individual_dedup_cnt = SEPARATE_BED.out.counts

        // Extract reads from merged qpass BAM matching the dedup BED → final BAM.
        ch_extract_in = RFC_DEDUP.out.bed
            .map { smeta, bed -> [smeta.id, bed] }
            .join(ch_merged_qpass_bam.map { smeta, bam, bai -> [smeta.id, bam] })
            .map { id, bed, bam -> [[id: id, strand: 'F'], bed, bam] }
        RFC_EXTRACT_DEDUP_READS(ch_extract_in)
        ch_final_bam        = RFC_EXTRACT_DEDUP_READS.out.bam
        ch_merged_dedup_cnt = RFC_EXTRACT_DEDUP_READS.out.counts
    }
    else if (dedup == 'umicollapse') {
        UMICOLLAPSE_DEDUP(ch_merged_qpass_bam)
        ch_final_bam        = UMICOLLAPSE_DEDUP.out.bam
        ch_merged_dedup_cnt = UMICOLLAPSE_DEDUP.out.counts

        // Split merged dedup BAM back to lanes → per-lane bam/bed/counts.
        ch_split_in = ch_lane_meta
            .combine(UMICOLLAPSE_DEDUP.out.bam.map { smeta, bam, bai -> [smeta.id, bam, bai] }, by: 0)
            .map { id, meta, bam, bai -> [meta, bam, bai] }
        SPLIT_DEDUP_BAM(ch_split_in)
        ch_individual_dedup_cnt = SPLIT_DEDUP_BAM.out.counts

        // Publish merged post-dedup BED (concat of per-lane post-dedup BEDs).
        CONCAT_POST_DEDUP_BED(
            SPLIT_DEDUP_BAM.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
    }

    DEEPTOOLS_BAMCOVERAGE(ch_final_bam)
    SPLIT_STRANDED_BAM(ch_final_bam)

    emit:
    transcriptome_bam       = STAR_ALIGN.out.transcriptome_bam // [ meta, txbam ] (only if emitted)
    genome_log              = STAR_ALIGN.out.log              // [ meta, log ]
    secondary_count         = STAR_ALIGN.out.secondary_count  // [ meta, count ]
    qpass_counts            = ch_qpass_counts                 // [ meta, t, p, s ]
    individual_dedup_counts = ch_individual_dedup_cnt         // [ meta, t, p, s ]
    merged_dedup_counts     = ch_merged_dedup_cnt             // [ smeta, t, p, s ] (empty if none)
}
