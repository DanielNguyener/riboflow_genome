// Transcriptome alignment stats: per-lane rows → combined essential + per-sample
// sums → published transcriptome_stats.csv / transcriptome_individual_stats.csv.
// Parallel to alignment_stats.nf (genome). Gated by params.transcriptome.run.
//
// Simplified vs genome stats: only total counts (no primary/secondary) are passed
// through. Per-lane totals sum correctly for transcriptome, so rfc update-dedup-counts
// is not needed (ext.use_merged_counts = false on TX_STATS_SUM).

include { TX_STATS_INDIVIDUAL }                         from '../../modules/local/tx_stats_individual.nf'
include { STATS_COMBINE as TX_COMBINE_INDIVIDUAL }      from '../../modules/local/stats_combine.nf'
include { STATS_COMBINE as TX_COMBINE_MERGED }          from '../../modules/local/stats_combine.nf'
include { STATS_SUM     as TX_STATS_SUM }               from '../../modules/local/stats_sum.nf'
include { STATS_PUBLISH as TX_STATS_PUBLISH }           from '../../modules/local/stats_publish.nf'

workflow TRANSCRIPTOME_STATS {
    take:
    ch_clip_log                // [ meta, log ]
    ch_filter_log              // [ meta, log ]
    ch_bowtie2_log             // [ meta, log ]
    ch_qpass_total_count       // [ meta, total ]
    ch_individual_dedup_counts // [ meta, total ] — total only

    main:
    def placeholder = file("${projectDir}/assets/NO_FILE.gz")

    // Per-lane stats input: join all per-lane channels on meta.
    ch_stats_in = ch_bowtie2_log
        .join(ch_clip_log)
        .join(ch_filter_log)
        .join(ch_qpass_total_count)
        .join(ch_individual_dedup_counts)
    TX_STATS_INDIVIDUAL(ch_stats_in)

    // Combined per-lane essential CSV.
    TX_COMBINE_INDIVIDUAL(TX_STATS_INDIVIDUAL.out.csv.map { meta, csv -> csv }.collect())

    // Per-sample sums — no merged-count override needed (per-lane totals are accurate).
    ch_grouped = TX_STATS_INDIVIDUAL.out.csv
        .map { meta, csv -> [meta.id, csv] }
        .groupTuple()
        .map { id, csvs -> [[id: id, strand: 'F'], csvs] }

    ch_sum_in = ch_grouped.map { smeta, csvs -> [smeta, csvs, placeholder, placeholder, placeholder] }
    TX_STATS_SUM(ch_sum_in)

    TX_COMBINE_MERGED(TX_STATS_SUM.out.csv.map { meta, csv -> csv }.collect())

    TX_STATS_PUBLISH(TX_COMBINE_INDIVIDUAL.out.csv, TX_COMBINE_MERGED.out.csv)

    emit:
    individual = TX_STATS_PUBLISH.out.individual
    merged     = TX_STATS_PUBLISH.out.merged
}
