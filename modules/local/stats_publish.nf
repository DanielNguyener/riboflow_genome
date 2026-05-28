// Publish the final stats CSVs. Ports `publish_stats` (RiboFlow.groovy:1758-1775).
process STATS_PUBLISH {
    executor 'local'

    input:
    path(individual_stats)
    path(merged_stats)

    output:
    path('individual_stats.csv'), emit: individual
    path('stats.csv'),            emit: merged

    script:
    """
    cp ${individual_stats} individual_stats.csv
    cp ${merged_stats} stats.csv
    """

    stub:
    """
    touch individual_stats.csv stats.csv
    """
}
