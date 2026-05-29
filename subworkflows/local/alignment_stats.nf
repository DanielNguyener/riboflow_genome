// Alignment stats: per-lane rows → combined essential + per-sample sums →
// published stats.csv / individual_stats.csv. (RiboFlow.groovy:1486-1775)

include { STATS_INDIVIDUAL }                        from '../../modules/local/stats_individual.nf'
include { STATS_COMBINE as COMBINE_INDIVIDUAL }     from '../../modules/local/stats_combine.nf'
include { STATS_COMBINE as COMBINE_MERGED }         from '../../modules/local/stats_combine.nf'
include { STATS_SUM }                               from '../../modules/local/stats_sum.nf'
include { STATS_PUBLISH }                           from '../../modules/local/stats_publish.nf'

workflow ALIGNMENT_STATS {
    take:
    ch_clip_log                // [ meta, log ]
    ch_filter_log              // [ meta, log ]
    ch_genome_log              // [ meta, log ]
    ch_secondary_count         // [ meta, count ]
    ch_qpass_total_count       // [ meta, total ]
    ch_qpass_primary_count     // [ meta, primary ]
    ch_qpass_secondary_count   // [ meta, secondary ]
    ch_qpass_unique_count      // [ meta, unique ] (empty if count_unique=false)
    ch_individual_dedup_counts // [ meta, t, p, s, u ]
    ch_merged_dedup_counts     // [ smeta, t, p, s, u ] (empty if dedup none)

    main:
    def dedup = Utils.resolve_dedup_method(params)
    def placeholder = file("${projectDir}/assets/NO_FILE.gz")

    // Per-lane stats input: join everything on the per-lane meta.
    ch_stats_in = ch_clip_log
        .join(ch_filter_log)
        .join(ch_genome_log)
        .join(ch_secondary_count)
        .join(ch_qpass_total_count)
        .join(ch_qpass_primary_count)
        .join(ch_qpass_secondary_count)
        .join(ch_individual_dedup_counts)
        .join(ch_qpass_unique_count)
    STATS_INDIVIDUAL(ch_stats_in)

    // Combined per-lane essential CSV.
    COMBINE_INDIVIDUAL(STATS_INDIVIDUAL.out.csv.map { meta, csv -> csv }.collect())

    // Per-sample sums.
    ch_grouped = STATS_INDIVIDUAL.out.csv
        .map { meta, csv -> [meta.id, csv] }
        .groupTuple()
        .map { id, csvs -> [[id: id, strand: 'F'], csvs] }

    if (dedup == 'position' || dedup == 'umicollapse') {
        ch_sum_in = ch_grouped
            .map { smeta, csvs -> [smeta.id, smeta, csvs] }
            .join(ch_merged_dedup_counts.map { smeta, t, p, s, u -> [smeta.id, t, p, s, u] })
            .map { id, smeta, csvs, t, p, s, u -> [smeta, csvs, t, p, s, u] }
    } else {
        ch_sum_in = ch_grouped.map { smeta, csvs -> [smeta, csvs, placeholder, placeholder, placeholder, placeholder] }
    }
    STATS_SUM(ch_sum_in)

    COMBINE_MERGED(STATS_SUM.out.csv.map { meta, csv -> csv }.collect())

    STATS_PUBLISH(COMBINE_INDIVIDUAL.out.csv, COMBINE_MERGED.out.csv)

    emit:
    individual = STATS_PUBLISH.out.individual
    merged     = STATS_PUBLISH.out.merged
}
