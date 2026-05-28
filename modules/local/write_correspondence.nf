// Publish the sample/lane → FASTQ correspondence table. Ports
// `write_fastq_correspondence` (RiboFlow.groovy:182-195).
process WRITE_CORRESPONDENCE {
    executor 'local'

    input:
    path(correspondence)

    output:
    path('index_fastq_correspondence.txt'), emit: txt

    script:
    """
    cat ${correspondence} > index_fastq_correspondence.txt
    """

    stub:
    """
    touch index_fastq_correspondence.txt
    """
}
