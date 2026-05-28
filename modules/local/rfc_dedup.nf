// Coordinate-based dedup via `rfc dedup`. Ports `genome_deduplicate_position`
// (RiboFlow.groovy:1069-1088) and `transcriptome_deduplicate_position` (:764-784).
process RFC_DEDUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${prefix}.post_dedup.bed"), emit: bed

    script:
    prefix = task.ext.prefix ?: "${meta.id}.genome"
    """
    sort -k1,1 -k2,2n -k3,3n ${bed} > sorted.bed
    rfc dedup -i sorted.bed -o ${prefix}.post_dedup.bed
    rm -f sorted.bed
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.genome"
    """
    touch ${prefix}.post_dedup.bed
    """
}
