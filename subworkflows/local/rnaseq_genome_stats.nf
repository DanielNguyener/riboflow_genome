// RNA-seq genome alignment stats. Parallel to alignment_stats.nf (ribo-seq):
// per-lane rows → combined essential + per-sample sums → published stats CSVs.
// Uses the RNA-seq dedup method and the rnaseq-specific stats modules.

include { RNASEQ_GENOME_STATS_INDIVIDUAL }                  from '../../modules/local/rnaseq_genome_stats_individual.nf'
include { STATS_COMBINE as RNASEQ_GENOME_COMBINE_INDIVIDUAL } from '../../modules/local/stats_combine.nf'
include { STATS_COMBINE as RNASEQ_GENOME_COMBINE_MERGED }     from '../../modules/local/stats_combine.nf'
include { STATS_SUM     as RNASEQ_GENOME_STATS_SUM }          from '../../modules/local/stats_sum.nf'
include { STATS_PUBLISH as RNASEQ_GENOME_STATS_PUBLISH }      from '../../modules/local/stats_publish.nf'

workflow RNASEQ_GENOME_STATS {
    take:
    ch_clip_log                // [ meta, log ]
    ch_filter_log              // [ meta, log ]
    ch_genome_log              // [ meta, log ]
    ch_secondary_count         // [ meta, count ]
    ch_qpass_total_count       // [ meta, total ]
    ch_qpass_primary_count     // [ meta, primary ]
    ch_qpass_secondary_count   // [ meta, secondary ]
    ch_qpass_unique_count      // [ meta, unique ]
    ch_individual_dedup_counts // [ meta, t, p, s, u ]
    ch_merged_dedup_counts     // [ smeta, t, p, s, u ] (empty if dedup none)

    main:
    def dedup = Utils.resolve_rnaseq_dedup_method(params)
    def placeholder = file("${projectDir}/assets/NO_FILE.gz")

    // When MAPQ > 0 (unique_only), SAMTOOLS_QPASS skips primary/secondary/unique
    // counts — those channels are empty and would break the join. Synthesise proxy
    // channels: primary = total (flag 2308 drops all secondaries so primary = total),
    // secondary = 0, unique = 0 (unused in unique_only mode but join needs a value).
    def mapq        = ((params.rnaseq?.genome?.mapping_quality_cutoff ?: params.rnaseq?.mapping_quality_cutoff ?: 4) as int)
    def unique_only = mapq > 0
    def zero_file   = file("${projectDir}/assets/zero.count")

    ch_qpass_primary_eff   = unique_only ? ch_qpass_total_count.map { m, f -> [m, f] }          : ch_qpass_primary_count
    ch_qpass_secondary_eff = unique_only ? ch_qpass_total_count.map { m, f -> [m, zero_file] }  : ch_qpass_secondary_count
    ch_qpass_unique_eff    = unique_only ? ch_qpass_total_count.map { m, f -> [m, zero_file] }  : ch_qpass_unique_count

    ch_stats_in = ch_clip_log
        .join(ch_filter_log)
        .join(ch_genome_log)
        .join(ch_secondary_count)
        .join(ch_qpass_total_count)
        .join(ch_qpass_primary_eff)
        .join(ch_qpass_secondary_eff)
        .join(ch_individual_dedup_counts)
        .join(ch_qpass_unique_eff)
    RNASEQ_GENOME_STATS_INDIVIDUAL(ch_stats_in)

    RNASEQ_GENOME_COMBINE_INDIVIDUAL(RNASEQ_GENOME_STATS_INDIVIDUAL.out.csv.map { meta, csv -> csv }.collect())

    ch_grouped = RNASEQ_GENOME_STATS_INDIVIDUAL.out.csv
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
    RNASEQ_GENOME_STATS_SUM(ch_sum_in)

    RNASEQ_GENOME_COMBINE_MERGED(RNASEQ_GENOME_STATS_SUM.out.csv.map { meta, csv -> csv }.collect())

    RNASEQ_GENOME_STATS_PUBLISH(RNASEQ_GENOME_COMBINE_INDIVIDUAL.out.csv, RNASEQ_GENOME_COMBINE_MERGED.out.csv)

    emit:
    individual = RNASEQ_GENOME_STATS_PUBLISH.out.individual
    merged     = RNASEQ_GENOME_STATS_PUBLISH.out.merged
}
