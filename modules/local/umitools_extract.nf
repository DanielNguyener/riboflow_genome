// UMI extraction (umicollapse path only). Ports `extract_umi_via_umi_tools`
// (RiboFlow.groovy:280-298).
process UMITOOLS_EXTRACT {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.umi_extracted.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.${meta.lane}.umi_extracted.log"),      emit: log

    script:
    def prefix = "${meta.id}.${meta.lane}"
    def args   = task.ext.args ?: (params.umi_tools_extract_arguments ?: '')
    """
    umi_tools extract -I ${fastq} -S ${prefix}.umi_extracted.fastq.gz \\
        -L ${prefix}.umi_extracted.log \\
        ${args}
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    echo | gzip -c > ${prefix}.umi_extracted.fastq.gz
    touch ${prefix}.umi_extracted.log
    """
}
