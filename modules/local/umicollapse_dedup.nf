// UMI-aware dedup. Ports `genome_deduplicate_umicollapse` (RiboFlow.groovy:1208-1246)
// and `transcriptome_deduplicate_umicollapse` (:916-941, no counts).
//
// ENV CONSOLIDATION: the original called a hand-shipped jar via a `java11`
// wrapper (`java11 -Xms512m -Xmx32g -Xss256m umicollapse.main.Main bam ...`).
// We now use the bioconda `umicollapse` CLI from the single conda env. Its
// wrapper only forwards `-Xm*` JVM flags as args and would mis-parse `-Xss`, so
// ALL JVM options are passed via `_JAVA_OPTIONS`; the CLI flags are unchanged.
process UMICOLLAPSE_DEDUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.merged_dedup.total.count"),
                     path("${meta.id}.merged_dedup.primary.count"),
                     path("${meta.id}.merged_dedup.secondary.count"),
                     path("${meta.id}.merged_dedup.unique.count"),
                     optional: true, emit: counts

    script:
    prefix          = task.ext.prefix ?: "${meta.id}.dedup"
    def args        = task.ext.args ?: (params.umicollapse_arguments ?: '')
    def jvm_opts    = task.ext.jvm_opts ?: '-Xms512m -Xmx32g -Xss256m'
    def emit_counts = task.ext.emit_counts ?: false
    def counts_cmd  = emit_counts ? """
    samtools view -@ ${task.cpus} -c        ${prefix}.bam > ${meta.id}.merged_dedup.total.count
    samtools view -@ ${task.cpus} -c -F 2304 ${prefix}.bam > ${meta.id}.merged_dedup.primary.count
    samtools view -@ ${task.cpus} -c -f 256  ${prefix}.bam > ${meta.id}.merged_dedup.secondary.count
    samtools view -@ ${task.cpus} -c -q 255  ${prefix}.bam > ${meta.id}.merged_dedup.unique.count
    """ : ''
    """
    _JAVA_OPTIONS="${jvm_opts}" umicollapse bam \\
        -i ${bam} \\
        -o ${prefix}.bam \\
        --umi-sep "_" \\
        --algo dir \\
        --merge mapqual \\
        --two-pass \\
        ${args}
    samtools index -@ ${task.cpus} ${prefix}.bam
    ${counts_cmd}
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}.dedup"
    def emit_counts = task.ext.emit_counts ?: false
    def counts_cmd  = emit_counts ? "echo 0 > ${meta.id}.merged_dedup.total.count; echo 0 > ${meta.id}.merged_dedup.primary.count; echo 0 > ${meta.id}.merged_dedup.secondary.count; echo 0 > ${meta.id}.merged_dedup.unique.count" : ''
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    ${counts_cmd}
    """
}
