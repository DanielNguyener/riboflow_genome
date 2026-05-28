// Extract reads from a qpass merged BAM matching an rfc-dedup BED →
// sample-level post-dedup BAM. With ext.emit_counts=true also emits merged
// total/primary/secondary counts. Ports
// `genome_convert_dedup_bed_to_bam_position` (RiboFlow.groovy:1162-1198) and
// `transcriptome_convert_dedup_bed_to_bam_position` (:859-883, no counts).
process RFC_EXTRACT_DEDUP_READS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(dedup_bed), path(qpass_bam)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.merged_dedup.total.count"),
                     path("${meta.id}.merged_dedup.primary.count"),
                     path("${meta.id}.merged_dedup.secondary.count"),
                     optional: true, emit: counts

    script:
    prefix          = task.ext.prefix ?: "${meta.id}.post_dedup"
    def emit_counts = task.ext.emit_counts ?: false
    def counts_cmd  = emit_counts ? """
    samtools view -@ ${task.cpus} -c        ${prefix}.bam > ${meta.id}.merged_dedup.total.count
    samtools view -@ ${task.cpus} -c -F 2304 ${prefix}.bam > ${meta.id}.merged_dedup.primary.count
    samtools view -@ ${task.cpus} -c -f 256  ${prefix}.bam > ${meta.id}.merged_dedup.secondary.count
    """ : ''
    """
    rfc extract-dedup-reads \\
        --bam ${qpass_bam} \\
        --bed ${dedup_bed} \\
        --output ${prefix}.bam
    samtools index -@ ${task.cpus} ${prefix}.bam
    ${counts_cmd}
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}.post_dedup"
    def emit_counts = task.ext.emit_counts ?: false
    def counts_cmd  = emit_counts ? "echo 0 > ${meta.id}.merged_dedup.total.count; echo 0 > ${meta.id}.merged_dedup.primary.count; echo 0 > ${meta.id}.merged_dedup.secondary.count" : ''
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    ${counts_cmd}
    """
}
