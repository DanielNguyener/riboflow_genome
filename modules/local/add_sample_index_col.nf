// Tag each per-lane BED record with the sample.lane id (column 7) so the merged
// rfc-dedup output can be split back per lane. Ports
// `add_sample_index_col_to_genome_bed` (RiboFlow.groovy:1012-1030).
process ADD_SAMPLE_INDEX_COL {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.genome.with_sample_index.bed"), emit: bed

    script:
    def prefix = "${meta.id}.${meta.lane}"
    """
    awk -v this_sample=${prefix} \\
        '{ print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"this_sample) }' ${bed} \\
        > ${prefix}.genome.with_sample_index.bed
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    touch ${prefix}.genome.with_sample_index.bed
    """
}
