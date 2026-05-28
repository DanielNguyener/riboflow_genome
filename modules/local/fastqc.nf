// Parametrized FastQC, gated on params.do_fastqc. One module replaces
// raw_fastqc / clipped_fastqc / genome_aligned_fastqc / genome_unaligned_fastqc
// (RiboFlow.groovy:233,312,494,523). `stage` (via ext.prefix) sets the output
// basename; the >20-byte guard skips empty FASTQs (matches the genome variants).
process FASTQC {
    tag "${prefix}"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*_fastqc.html"), path("*_fastqc.zip"), emit: report

    when:
    params.do_fastqc

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.lane}"
    """
    ln -sf ${fastq} ${prefix}.fastq.gz
    if [ \$(stat -L -c%s ${prefix}.fastq.gz) -gt 20 ]; then
        fastqc ${prefix}.fastq.gz --outdir=\$PWD -t ${task.cpus}
    else
        echo "File is empty, skipping FastQC"
        touch ${prefix}_fastqc.html ${prefix}_fastqc.zip
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.lane}"
    """
    touch ${prefix}_fastqc.html ${prefix}_fastqc.zip
    """
}
