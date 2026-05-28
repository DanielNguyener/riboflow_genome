// Combine per-row stats CSVs and add percentage rows. Ports
// `combine_individual_genome_alignment_stats` (RiboFlow.groovy:1627-1654) and
// `combine_merged_genome_alignment_stats` (:1718-1744). ext.prefix sets the
// output basename (genome_individual_essential / genome_merged_essential).
process STATS_COMBINE {
    executor 'local'

    input:
    path(stat_tables)

    output:
    path("${prefix}.csv"), emit: csv

    script:
    prefix = task.ext.prefix ?: 'genome_individual_essential'
    def n  = stat_tables instanceof List ? stat_tables.size() : 1
    if (n == 0) {
        """
        echo "No statistics data available" > ${prefix}.csv
        """
    } else {
        """
        rfc merge overall-stats -o raw_${prefix}.csv ${stat_tables}
        rfc genome-stats-percentage -i raw_${prefix}.csv -o ${prefix}.csv
        """
    }

    stub:
    prefix = task.ext.prefix ?: 'genome_individual_essential'
    """
    touch ${prefix}.csv
    """
}
