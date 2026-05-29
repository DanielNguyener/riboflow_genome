// Quality-pass filter (MAPQ + FLAG mask) with per-type count files.
// Ports `genome_quality_filter` (RiboFlow.groovy:576-609) and, with
// ext.presort=true, `transcriptome_sort_and_filter` (:675-701) — the STAR
// transcriptome BAM is QNAME-sorted so it must be coord-sorted after filtering.
//
// ext.emit_primary_secondary (default true): when false, skip primary/secondary
// count commands (transcriptome path only needs total).
// ext.count_unique (default false): emit MAPQ-255 unique count (genome multi mode).
process SAMTOOLS_QPASS {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.qpass.bam"), path("${prefix}.qpass.bam.bai"), emit: bam
    tuple val(meta), path("${prefix}.qpass.total.count"),   emit: total_count
    tuple val(meta), path("${prefix}.qpass.primary.count"), optional: true, emit: primary_count
    tuple val(meta), path("${prefix}.qpass.secondary.count"), optional: true, emit: secondary_count
    tuple val(meta), path("${prefix}.qpass.unique.count"),  optional: true, emit: unique_count

    script:
    prefix                  = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def presort             = task.ext.presort ?: false
    def mapq                = (task.ext.mapq != null) ? task.ext.mapq : params.genome.mapping_quality_cutoff
    def emit_primary_second = (task.ext.emit_primary_secondary != null) ? task.ext.emit_primary_secondary : true
    def count_unique        = task.ext.count_unique ?: false
    def filter_flags        = (task.ext.filter_flags != null) ? task.ext.filter_flags : params.genome.ribo_filter_flags
    def sort_threads        = Math.min(task.cpus as int, 8)
    def sort_mem            = Utils.samtools_sort_mem_per_thread_mb(task)
    def make_bam            = presort \
        ? "samtools view -h -bq ${mapq} -F ${filter_flags} ${bam} | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${prefix}.qpass.bam -" \
        : "samtools view -@ ${task.cpus} -bq ${mapq} -F ${filter_flags} ${bam} > ${prefix}.qpass.bam"
    def ps_cmd  = emit_primary_second ? """
    samtools view -@ ${task.cpus} -c -F 2304 ${prefix}.qpass.bam > ${prefix}.qpass.primary.count
    samtools view -@ ${task.cpus} -c -f 256  ${prefix}.qpass.bam > ${prefix}.qpass.secondary.count
    """ : ''
    def uniq_cmd = count_unique ? "samtools view -@ ${task.cpus} -c -q 255 ${prefix}.qpass.bam > ${prefix}.qpass.unique.count" : ''
    """
    ${make_bam}
    samtools index -@ ${task.cpus} ${prefix}.qpass.bam
    samtools view -@ ${task.cpus} -c ${prefix}.qpass.bam > ${prefix}.qpass.total.count
    ${ps_cmd}
    ${uniq_cmd}
    """

    stub:
    prefix                  = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def emit_primary_second = (task.ext.emit_primary_secondary != null) ? task.ext.emit_primary_secondary : true
    def count_unique        = task.ext.count_unique ?: false
    def ps_cmd   = emit_primary_second ? "echo 0 > ${prefix}.qpass.primary.count; echo 0 > ${prefix}.qpass.secondary.count" : ''
    def uniq_cmd = count_unique ? "echo 0 > ${prefix}.qpass.unique.count" : ''
    """
    touch ${prefix}.qpass.bam ${prefix}.qpass.bam.bai
    echo 0 > ${prefix}.qpass.total.count
    ${ps_cmd}
    ${uniq_cmd}
    """
}
