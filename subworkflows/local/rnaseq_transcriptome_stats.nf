// RNA-seq transcriptome alignment stats. Parallel to transcriptome_stats.nf
// (ribo-seq). Only total counts pass through; per-lane totals sum correctly so no
// merged-count override is needed (use_merged_counts=false on RNASEQ_TX_STATS_SUM).

include { RNASEQ_TX_STATS_INDIVIDUAL }                  from '../../modules/local/rnaseq_tx_stats_individual.nf'
include { STATS_COMBINE as RNASEQ_TX_COMBINE_INDIVIDUAL } from '../../modules/local/stats_combine.nf'
include { STATS_COMBINE as RNASEQ_TX_COMBINE_MERGED }     from '../../modules/local/stats_combine.nf'
include { STATS_SUM     as RNASEQ_TX_STATS_SUM }          from '../../modules/local/stats_sum.nf'
include { STATS_PUBLISH as RNASEQ_TX_STATS_PUBLISH }      from '../../modules/local/stats_publish.nf'

workflow RNASEQ_TX_STATS {
    take:
    ch_clip_log                // [ meta, log ]
    ch_filter_log              // [ meta, log ]
    ch_bowtie2_log             // [ meta, log ]
    ch_qpass_total_count       // [ meta, total ]
    ch_individual_dedup_counts // [ meta, total ] — total only

    main:
    def placeholder = file("${projectDir}/assets/NO_FILE.gz")

    ch_stats_in = ch_bowtie2_log
        .join(ch_clip_log)
        .join(ch_filter_log)
        .join(ch_qpass_total_count)
        .join(ch_individual_dedup_counts)
    RNASEQ_TX_STATS_INDIVIDUAL(ch_stats_in)

    RNASEQ_TX_COMBINE_INDIVIDUAL(RNASEQ_TX_STATS_INDIVIDUAL.out.csv.map { meta, csv -> csv }.collect())

    ch_grouped = RNASEQ_TX_STATS_INDIVIDUAL.out.csv
        .map { meta, csv -> [meta.id, csv] }
        .groupTuple()
        .map { id, csvs -> [[id: id, strand: 'F'], csvs] }

    ch_sum_in = ch_grouped.map { smeta, csvs -> [smeta, csvs, placeholder, placeholder, placeholder, placeholder] }
    RNASEQ_TX_STATS_SUM(ch_sum_in)

    RNASEQ_TX_COMBINE_MERGED(RNASEQ_TX_STATS_SUM.out.csv.map { meta, csv -> csv }.collect())

    RNASEQ_TX_STATS_PUBLISH(RNASEQ_TX_COMBINE_INDIVIDUAL.out.csv, RNASEQ_TX_COMBINE_MERGED.out.csv)

    emit:
    individual = RNASEQ_TX_STATS_PUBLISH.out.individual
    merged     = RNASEQ_TX_STATS_PUBLISH.out.merged
}
