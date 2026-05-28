// Quality-pass filter (MAPQ + FLAG mask) with total/primary/secondary counts.
// Ports `genome_quality_filter` (RiboFlow.groovy:576-609) and, with
// ext.presort=true, `transcriptome_sort_and_filter` (:675-701) — the STAR
// transcriptome BAM is QNAME-sorted so it must be coord-sorted after filtering.
process SAMTOOLS_QPASS {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.qpass.bam"), path("${prefix}.qpass.bam.bai"), emit: bam
    tuple val(meta), path("${prefix}.qpass.total.count"),
                     path("${prefix}.qpass.primary.count"),
                     path("${prefix}.qpass.secondary.count"), emit: counts
    tuple val(meta), path("${prefix}.qpass.unique.count"), optional: true, emit: unique_count

    script:
    prefix              = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def presort         = task.ext.presort ?: false
    def mapq            = (task.ext.mapq != null) ? task.ext.mapq : params.mapping_quality_cutoff
    def count_unique    = task.ext.count_unique ?: false
    def sort_threads    = Math.min(task.cpus as int, 8)
    def sort_mem        = Utils.samtools_sort_mem_per_thread_mb(task)
    def make_bam        = presort \
        ? "samtools view -h -bq ${mapq} -F ${params.ribo_filter_flags} ${bam} | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${prefix}.qpass.bam -" \
        : "samtools view -@ ${task.cpus} -bq ${mapq} -F ${params.ribo_filter_flags} ${bam} > ${prefix}.qpass.bam"
    def unique_cmd      = count_unique ? "samtools view -@ ${task.cpus} -c -q 255 ${prefix}.qpass.bam > ${prefix}.qpass.unique.count" : ''
    """
    ${make_bam}
    samtools index -@ ${task.cpus} ${prefix}.qpass.bam
    samtools view -@ ${task.cpus} -c        ${prefix}.qpass.bam > ${prefix}.qpass.total.count
    samtools view -@ ${task.cpus} -c -F 2304 ${prefix}.qpass.bam > ${prefix}.qpass.primary.count
    samtools view -@ ${task.cpus} -c -f 256  ${prefix}.qpass.bam > ${prefix}.qpass.secondary.count
    ${unique_cmd}
    """

    stub:
    prefix           = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def count_unique = task.ext.count_unique ?: false
    def unique_cmd   = count_unique ? "echo 0 > ${prefix}.qpass.unique.count" : ''
    """
    touch ${prefix}.qpass.bam ${prefix}.qpass.bam.bai
    echo 0 > ${prefix}.qpass.total.count
    echo 0 > ${prefix}.qpass.primary.count
    echo 0 > ${prefix}.qpass.secondary.count
    ${unique_cmd}
    """
}
