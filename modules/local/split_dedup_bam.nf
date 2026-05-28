// Split a merged dedup BAM back to one lane by read group. With
// ext.emit_bed_counts=true also emits the per-lane BED and total/primary/
// secondary counts. Ports `split_genome_dedup_bam_to_individual`
// (RiboFlow.groovy:1264-1309) and `split_transcriptome_dedup_bam_to_individual`
// (:952-975, bam only).
process SPLIT_DEDUP_BAM {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(merged_bam), path(merged_bai)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.genome.post_dedup.bed"), optional: true, emit: bed
    tuple val(meta), path("${meta.id}.${meta.lane}.dedup.total.count"),
                     path("${meta.id}.${meta.lane}.dedup.primary.count"),
                     path("${meta.id}.${meta.lane}.dedup.secondary.count"),
                     optional: true, emit: counts

    script:
    prefix          = task.ext.prefix ?: "${meta.id}.${meta.lane}.post_dedup"
    def s           = "${meta.id}.${meta.lane}"
    def emit_extra  = task.ext.emit_bed_counts ?: false
    def extra_cmd   = emit_extra ? """
    samtools view -@ ${task.cpus} -c        ${prefix}.bam > ${s}.dedup.total.count
    samtools view -@ ${task.cpus} -c -F 2304 ${prefix}.bam > ${s}.dedup.primary.count
    samtools view -@ ${task.cpus} -c -f 256  ${prefix}.bam > ${s}.dedup.secondary.count
    if [ \$(cat ${s}.dedup.total.count) -eq 0 ]; then
        touch ${s}.genome.post_dedup.bed
    else
        bamToBed -i ${prefix}.bam > ${s}.genome.post_dedup.bed
    fi
    """ : ''
    """
    if ! samtools view -H ${merged_bam} | grep -q "^@RG"; then
        echo "ERROR: No read groups found in ${merged_bam}" >&2; exit 1
    fi
    samtools view -@ ${task.cpus} -B -r ${s} ${merged_bam} -o ${prefix}.bam
    samtools index -@ ${task.cpus} ${prefix}.bam
    ${extra_cmd}
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}.${meta.lane}.post_dedup"
    def s           = "${meta.id}.${meta.lane}"
    def emit_extra  = task.ext.emit_bed_counts ?: false
    def extra_cmd   = emit_extra ? "touch ${s}.genome.post_dedup.bed; echo 0 > ${s}.dedup.total.count; echo 0 > ${s}.dedup.primary.count; echo 0 > ${s}.dedup.secondary.count" : ''
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    ${extra_cmd}
    """
}
