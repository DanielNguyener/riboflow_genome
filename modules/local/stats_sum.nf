// Sum per-lane stats into a per-sample row. Ports
// `sum_individual_genome_alignment_stats` (RiboFlow.groovy:1677-1710). When
// ext.use_merged_counts=true (dedup position/umicollapse) the merged dedup
// counts override the summed per-lane dedup rows via `rfc update-dedup-counts`;
// otherwise (dedup none) the placeholder count files are ignored.
process STATS_SUM {
    executor 'local'
    tag "${meta.id}"

    input:
    tuple val(meta), path(stat_files),
          path(dedup_total,     stageAs: 'merged_dedup.total.count'),
          path(dedup_primary,   stageAs: 'merged_dedup.primary.count'),
          path(dedup_secondary, stageAs: 'merged_dedup.secondary.count')

    output:
    tuple val(meta), path("${meta.id}.${task.ext.label ?: 'genome'}_merged.csv"), emit: csv

    script:
    def label      = task.ext.label ?: 'genome'
    def use_counts = task.ext.use_merged_counts ?: false
    if (use_counts) {
        """
        rfc sum-stats -n ${meta.id} -o ${meta.id}.${label}_merged.tmp.csv ${stat_files}
        rfc update-dedup-counts \\
            --dedup-total-file merged_dedup.total.count \\
            --dedup-primary-file merged_dedup.primary.count \\
            --dedup-secondary-file merged_dedup.secondary.count \\
            --input-csv ${meta.id}.${label}_merged.tmp.csv \\
            --output-csv ${meta.id}.${label}_merged.csv
        """
    } else {
        """
        rfc sum-stats -n ${meta.id} -o ${meta.id}.${label}_merged.csv ${stat_files}
        """
    }

    stub:
    def label = task.ext.label ?: 'genome'
    """
    printf ',${meta.id}\\ntotal_reads,0\\n' > ${meta.id}.${label}_merged.csv
    """
}
