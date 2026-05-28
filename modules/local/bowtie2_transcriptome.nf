// Bowtie2 alignment against transcriptome index. Ports `transcriptome_alignment`
// from upstream RiboFlow (ribosomeprofiling/riboflow RiboFlow.groovy).
// Outputs a coord-sorted BAM + bowtie2 log. No aligned/unaligned FASTQ —
// the aligned BAM is what feeds the transcriptome dedup/ribopy chain.
process BOWTIE2_TRANSCRIPTOME {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)
    tuple val(index_base), path(index_files)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.bam"), emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.log"), emit: log

    script:
    def prefix       = "${meta.id}.${meta.lane}"
    def args         = task.ext.args ?: (params.alignment_arguments?.transcriptome ?: '')
    def aln_threads  = Math.min(task.cpus as int, 16)
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = Utils.samtools_sort_mem_per_thread_mb(task)
    """
    set -o pipefail
    bowtie2 ${args} \\
            -x ${index_base} -q ${fastq} \\
            --threads ${aln_threads} \\
            --rg-id "${prefix}" \\
            --rg "SM:${meta.id}" \\
            --rg "LB:${prefix}" \\
            --rg "PL:ILLUMINA" \\
                     2> ${prefix}.transcriptome_alignment.log \\
            | samtools sort -@ ${sort_threads} -m ${sort_mem}M \\
                     -o ${prefix}.transcriptome_alignment.bam -
    samtools index -@ ${sort_threads} ${prefix}.transcriptome_alignment.bam
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    touch ${prefix}.transcriptome_alignment.bam
    touch ${prefix}.transcriptome_alignment.log
    """
}
