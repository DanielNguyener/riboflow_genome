// Split a merged post-dedup BED back to a single lane by the sample.lane id in
// column 7. With ext.emit_counts=true also emits per-lane total/primary/secondary
// counts. Ports `separate_genome_bed_post_dedup` (RiboFlow.groovy:1120-1151) and
// `separate_transcriptome_bed_post_dedup` (:801-819, no counts).
process SEPARATE_BED {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(merged_bed)

    output:
    tuple val(meta), path("${prefix}.post_dedup.bed"), emit: bed
    tuple val(meta), path("${meta.id}.${meta.lane}.dedup.total.count"),
                     path("${meta.id}.${meta.lane}.dedup.primary.count"),
                     path("${meta.id}.${meta.lane}.dedup.secondary.count"),
                     path("${meta.id}.${meta.lane}.dedup.unique.count"),
                     optional: true, emit: counts

    script:
    prefix          = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome"
    def s           = "${meta.id}.${meta.lane}"
    def emit_counts = task.ext.emit_counts ?: false
    def counts_cmd  = emit_counts ? """
    total=\$(wc -l < ${prefix}.post_dedup.bed)
    primary=\$(awk '{print \$4}' ${prefix}.post_dedup.bed | sort -u | wc -l)
    secondary=\$((total - primary))
    unique=\$(awk '\$5 >= 255' ${prefix}.post_dedup.bed | wc -l)
    echo \${total}     > ${s}.dedup.total.count
    echo \${primary}   > ${s}.dedup.primary.count
    echo \${secondary} > ${s}.dedup.secondary.count
    echo \${unique}    > ${s}.dedup.unique.count
    """ : ''
    """
    awk -v s=${s} \\
        '\$7 == s { print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6 }' \\
        ${merged_bed} > ${prefix}.post_dedup.bed
    ${counts_cmd}
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome"
    def s           = "${meta.id}.${meta.lane}"
    def emit_counts = task.ext.emit_counts ?: false
    def counts_cmd  = emit_counts ? "echo 0 > ${s}.dedup.total.count; echo 0 > ${s}.dedup.primary.count; echo 0 > ${s}.dedup.secondary.count; echo 0 > ${s}.dedup.unique.count" : ''
    """
    touch ${prefix}.post_dedup.bed
    ${counts_cmd}
    """
}
