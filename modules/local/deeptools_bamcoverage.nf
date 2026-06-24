// Strand-specific bigWigs. Ports `genome_create_strand_specific_bigwigs`
// (RiboFlow.groovy:1391-1430). deepTools --filterRNAstrand assumes reverse-
// stranded libraries, so for forward-stranded ribo-seq (strand F/FR) the
// plus/minus filters are swapped.
process DEEPTOOLS_BAMCOVERAGE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    val assay   // filename label: 'ribo' for ribo-seq, 'rna' for RNA-seq

    output:
    tuple val(meta), path("${meta.id}.${assay}.plus.bigWig"),
                     path("${meta.id}.${assay}.minus.bigWig"), emit: bigwig

    when:
    params.get('do_bigwig', false)

    script:
    def strand_arg = meta.strand ?: 'F'
    def bw_threads = Math.min(task.cpus as int, 8)
    """
    if [ "${strand_arg}" == "F" ] || [ "${strand_arg}" == "FR" ]; then
        bamCoverage -b ${bam} -o ${meta.id}.${assay}.plus.bigWig \\
            --filterRNAstrand reverse --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
        bamCoverage -b ${bam} -o ${meta.id}.${assay}.minus.bigWig \\
            --filterRNAstrand forward --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
    else
        bamCoverage -b ${bam} -o ${meta.id}.${assay}.plus.bigWig \\
            --filterRNAstrand forward --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
        bamCoverage -b ${bam} -o ${meta.id}.${assay}.minus.bigWig \\
            --filterRNAstrand reverse --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
    fi
    """

    stub:
    """
    touch ${meta.id}.${assay}.plus.bigWig ${meta.id}.${assay}.minus.bigWig
    """
}
