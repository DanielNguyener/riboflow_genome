// BAM → BED with an empty-input guard. Ports `individual_genome_bam_to_bed`
// (RiboFlow.groovy:983-1002) and, with ext.append_index=true,
// `transcriptome_qpass_bam_to_bed` (:715-734) which tags each record with the
// sample.lane id in column 7.
process BAM_TO_BED {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed

    script:
    prefix       = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome.qpass"
    def s        = "${meta.id}.${meta.lane}"
    def to_bed   = (task.ext.append_index ?: false) \
        ? "bamToBed -i ${bam} | awk -v s=${s} '{ print \$0\"\\t\"s }' > ${prefix}.bed" \
        : "bamToBed -i ${bam} > ${prefix}.bed"
    """
    if [ \$(samtools view -@ ${task.cpus} -c ${bam}) -eq 0 ]; then
        touch ${prefix}.bed
    else
        ${to_bed}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome.qpass"
    """
    touch ${prefix}.bed
    """
}
