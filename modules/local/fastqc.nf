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
    def reads = (fastq instanceof List) ? fastq : [fastq]
    if (reads.size() > 1) {
        // Paired-end: FastQC each mate separately, tagged _R1 / _R2.
        """
        idx=0
        for f in ${reads.join(' ')}; do
            idx=\$((idx+1))
            nm=${prefix}_R\${idx}
            ln -sf \$f \${nm}.fastq.gz
            if [ \$(stat -L -c%s \${nm}.fastq.gz) -gt 20 ]; then
                fastqc \${nm}.fastq.gz --outdir=\$PWD -t ${task.cpus}
            else
                echo "File is empty, skipping FastQC"
                touch \${nm}_fastqc.html \${nm}_fastqc.zip
            fi
        done
        """
    } else {
        """
        ln -sf ${fastq} ${prefix}.fastq.gz
        if [ \$(stat -L -c%s ${prefix}.fastq.gz) -gt 20 ]; then
            fastqc ${prefix}.fastq.gz --outdir=\$PWD -t ${task.cpus}
        else
            echo "File is empty, skipping FastQC"
            touch ${prefix}_fastqc.html ${prefix}_fastqc.zip
        fi
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.lane}"
    def reads = (fastq instanceof List) ? fastq : [fastq]
    if (reads.size() > 1) {
        """
        touch ${prefix}_R1_fastqc.html ${prefix}_R1_fastqc.zip ${prefix}_R2_fastqc.html ${prefix}_R2_fastqc.zip
        """
    } else {
        """
        touch ${prefix}_fastqc.html ${prefix}_fastqc.zip
        """
    }
}
