// Split the final sample BAM into +/- strand BAM/BED. Gated on do_strand_split.
// Ports `genome_split_stranded_bam` (RiboFlow.groovy:1432-1464).
process SPLIT_STRANDED_BAM {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.ribo.plus.bam"),  path("${meta.id}.ribo.plus.bam.bai"),
                     path("${meta.id}.ribo.minus.bam"), path("${meta.id}.ribo.minus.bam.bai"),
                     path("${meta.id}.ribo.plus.bed"),  path("${meta.id}.ribo.minus.bed"), emit: stranded

    when:
    params.get('do_strand_split', false)

    script:
    def strand_arg = meta.strand ?: 'F'
    """
    if [ "${strand_arg}" == "F" ] || [ "${strand_arg}" == "FR" ]; then
        samtools view -@ ${task.cpus} -b -F 2064        -o ${meta.id}.ribo.plus.bam  ${bam}
        samtools view -@ ${task.cpus} -b -f 16 -F 2048  -o ${meta.id}.ribo.minus.bam ${bam}
    else
        samtools view -@ ${task.cpus} -b -f 16 -F 2048  -o ${meta.id}.ribo.plus.bam  ${bam}
        samtools view -@ ${task.cpus} -b -F 2064        -o ${meta.id}.ribo.minus.bam ${bam}
    fi
    samtools index -@ ${task.cpus} ${meta.id}.ribo.plus.bam
    samtools index -@ ${task.cpus} ${meta.id}.ribo.minus.bam
    bamToBed -i ${meta.id}.ribo.plus.bam  > ${meta.id}.ribo.plus.bed
    bamToBed -i ${meta.id}.ribo.minus.bam > ${meta.id}.ribo.minus.bed
    """

    stub:
    """
    touch ${meta.id}.ribo.plus.bam ${meta.id}.ribo.plus.bam.bai
    touch ${meta.id}.ribo.minus.bam ${meta.id}.ribo.minus.bam.bai
    touch ${meta.id}.ribo.plus.bed ${meta.id}.ribo.minus.bed
    """
}
