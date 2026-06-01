// STAR genome alignment for the RNA-seq path. SE/PE aware. Mirrors star_align.nf
// but uses params.star.rnaseq_arguments, never emits a transcriptome BAM, and
// names outputs with the rnaseq.genome_alignment.* prefix (storeDir separation).
process STAR_ALIGN_RNASEQ {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(reads)
    path genome_dir

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.bam"),                 emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.aligned.fastq.gz"),    emit: aligned
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.unaligned.fastq.gz"),  emit: unaligned
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.log"),                 emit: log
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.genome_alignment.secondary.count"),     emit: secondary_count

    script:
    def prefix       = "${meta.id}.${meta.lane}"
    def is_pe        = meta.single_end == false
    def reads_in     = is_pe ? "${reads[0]} ${reads[1]}" : "${reads}"
    def star_args    = task.ext.star_args ?: params.star?.rnaseq_arguments
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = Utils.samtools_sort_mem_per_thread_mb(task)
    """
    set -o pipefail
    mkdir -p star_out
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${reads_in} \\
        --readFilesCommand zcat \\
        --readFilesType Fastx \\
        --genomeDir ${genome_dir} \\
        --genomeLoad NoSharedMemory \\
        ${star_args} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes All \\
        --outSAMstrandField intronMotif \\
        --outSAMattrRGline ID:${prefix} SM:${meta.id} PL:ILLUMINA \\
        --outReadsUnmapped Fastx \\
        --outFileNamePrefix star_out/

    samtools sort -@ ${sort_threads} -m ${sort_mem}M \\
        -o ${prefix}.rnaseq.genome_alignment.bam \\
        star_out/Aligned.sortedByCoord.out.bam
    rm -f star_out/Aligned.sortedByCoord.out.bam
    samtools index -@ ${sort_threads} ${prefix}.rnaseq.genome_alignment.bam

    samtools fastq -@ ${task.cpus} -F 4 ${prefix}.rnaseq.genome_alignment.bam \\
        | gzip > ${prefix}.rnaseq.genome_alignment.aligned.fastq.gz

    if [ -f star_out/Unmapped.out.mate1 ]; then
        gzip -c star_out/Unmapped.out.mate1 > ${prefix}.rnaseq.genome_alignment.unaligned.fastq.gz
    else
        echo -n | gzip > ${prefix}.rnaseq.genome_alignment.unaligned.fastq.gz
    fi

    cp star_out/Log.final.out ${prefix}.rnaseq.genome_alignment.log

    samtools view -@ ${task.cpus} -c -f 256 ${prefix}.rnaseq.genome_alignment.bam \\
        > ${prefix}.rnaseq.genome_alignment.secondary.count
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    touch ${prefix}.rnaseq.genome_alignment.bam
    echo | gzip -c > ${prefix}.rnaseq.genome_alignment.aligned.fastq.gz
    echo | gzip -c > ${prefix}.rnaseq.genome_alignment.unaligned.fastq.gz
    touch ${prefix}.rnaseq.genome_alignment.log
    echo 0 > ${prefix}.rnaseq.genome_alignment.secondary.count
    """
}
