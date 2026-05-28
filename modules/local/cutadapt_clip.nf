// Adapter trimming. Ports `clip` (RiboFlow.groovy:257-271).
process CUTADAPT_CLIP {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.clipped.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.${meta.lane}.clipped.log"),      emit: log

    script:
    def prefix = "${meta.id}.${meta.lane}"
    def args   = task.ext.args ?: params.clip_arguments
    """
    cutadapt --cores=${task.cpus} ${args} ${fastq} 2>${prefix}.clipped.log \\
        | gzip -c > ${prefix}.clipped.fastq.gz
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    echo | gzip -c > ${prefix}.clipped.fastq.gz
    touch ${prefix}.clipped.log
    """
}
