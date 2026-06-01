// STAR genome alignment. Ports `genome_alignment` (RiboFlow.groovy:401-483).
// Optionally emits a transcriptome-coordinate BAM (star.output_transcriptome_bam).
process STAR_ALIGN {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)
    path genome_dir

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.bam"),                   emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.bam"), optional: true, emit: transcriptome_bam
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.aligned.fastq.gz"),      emit: aligned
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.unaligned.fastq.gz"),    emit: unaligned
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.log"),                   emit: log
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_alignment.secondary.count"),       emit: secondary_count

    script:
    def prefix         = "${meta.id}.${meta.lane}"
    def emit_tx_bam    = (params.star ?: [:]).output_transcriptome_bam ?: false
    def quant_mode_arg = emit_tx_bam ? '--quantMode TranscriptomeSAM \\\n        ' : ''
    def tx_bam_cmd     = emit_tx_bam \
        ? "mv star_out/Aligned.toTranscriptome.out.bam ${prefix}.transcriptome_alignment.bam" \
        : ''
    def sort_threads   = Math.min(task.cpus as int, 8)
    def sort_mem       = Utils.samtools_sort_mem_per_thread_mb(task)
    def star_args      = task.ext.star_args ?: params.star.ribo_arguments
    """
    set -o pipefail
    mkdir -p star_out
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${fastq} \\
        --readFilesCommand zcat \\
        --readFilesType Fastx \\
        --genomeDir ${genome_dir} \\
        --genomeLoad NoSharedMemory \\
        ${star_args} \\
        --outSAMtype BAM SortedByCoordinate \\
        ${quant_mode_arg}--outSAMattributes All \\
        --outSAMstrandField intronMotif \\
        --outSAMattrRGline ID:${prefix} SM:${meta.id} PL:ILLUMINA \\
        --outReadsUnmapped Fastx \\
        --outFileNamePrefix star_out/

    samtools sort -@ ${sort_threads} -m ${sort_mem}M \\
        -o ${prefix}.genome_alignment.bam \\
        star_out/Aligned.sortedByCoord.out.bam
    rm -f star_out/Aligned.sortedByCoord.out.bam
    samtools index -@ ${sort_threads} ${prefix}.genome_alignment.bam

    ${tx_bam_cmd}

    samtools fastq -@ ${task.cpus} -F 4 ${prefix}.genome_alignment.bam \\
        | gzip > ${prefix}.genome_alignment.aligned.fastq.gz

    if [ -f star_out/Unmapped.out.mate1 ]; then
        gzip -c star_out/Unmapped.out.mate1 > ${prefix}.genome_alignment.unaligned.fastq.gz
    else
        echo -n | gzip > ${prefix}.genome_alignment.unaligned.fastq.gz
    fi

    cp star_out/Log.final.out ${prefix}.genome_alignment.log

    samtools view -@ ${task.cpus} -c -f 256 ${prefix}.genome_alignment.bam \\
        > ${prefix}.genome_alignment.secondary.count
    """

    stub:
    def prefix      = "${meta.id}.${meta.lane}"
    def emit_tx_bam = (params.star ?: [:]).output_transcriptome_bam ?: false
    def tx_stub     = emit_tx_bam ? "touch ${prefix}.transcriptome_alignment.bam" : ''
    """
    touch ${prefix}.genome_alignment.bam
    ${tx_stub}
    echo | gzip -c > ${prefix}.genome_alignment.aligned.fastq.gz
    echo | gzip -c > ${prefix}.genome_alignment.unaligned.fastq.gz
    touch ${prefix}.genome_alignment.log
    echo 0 > ${prefix}.genome_alignment.secondary.count
    """
}
