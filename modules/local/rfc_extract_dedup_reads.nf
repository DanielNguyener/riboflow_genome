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
                     optional: true, emit: total_count
    tuple val(meta), path("${meta.id}.merged_dedup.primary.count"),
                     path("${meta.id}.merged_dedup.secondary.count"),
                     path("${meta.id}.merged_dedup.unique.count"),
                     optional: true, emit: detail_counts

    script:
    prefix               = task.ext.prefix ?: "${meta.id}.post_dedup"
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : true
    def total_cmd  = emit_counts ? "samtools view -@ ${task.cpus} -c ${prefix}.bam > ${meta.id}.merged_dedup.total.count" : ''
    def detail_cmd = (emit_counts && emit_full_counts) ? """
    samtools view -@ ${task.cpus} -c -F 2304 ${prefix}.bam > ${meta.id}.merged_dedup.primary.count
    samtools view -@ ${task.cpus} -c -f 256  ${prefix}.bam > ${meta.id}.merged_dedup.secondary.count
    samtools view -@ ${task.cpus} -c -q 255  ${prefix}.bam > ${meta.id}.merged_dedup.unique.count
    """ : ''
    """
    if [ -s ${dedup_bed} ]; then
        rfc extract-dedup-reads \\
            --bam ${qpass_bam} \\
            --bed ${dedup_bed} \\
            --output ${prefix}.bam
    else
        samtools view -H ${qpass_bam} | samtools view -b -o ${prefix}.bam
    fi
    samtools index -@ ${task.cpus} ${prefix}.bam
    ${total_cmd}
    ${detail_cmd}
    """

    stub:
    prefix               = task.ext.prefix ?: "${meta.id}.post_dedup"
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : true
    def total_cmd  = emit_counts ? "echo 0 > ${meta.id}.merged_dedup.total.count" : ''
    def detail_cmd = (emit_counts && emit_full_counts) ? "echo 0 > ${meta.id}.merged_dedup.primary.count; echo 0 > ${meta.id}.merged_dedup.secondary.count; echo 0 > ${meta.id}.merged_dedup.unique.count" : ''
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    ${total_cmd}
    ${detail_cmd}
    """
}
