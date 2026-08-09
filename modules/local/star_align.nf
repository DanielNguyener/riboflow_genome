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
    def star_args      = task.ext.star_args ?: params.star.ribo_arguments
    // STAR holds the genome index (~40 GB) in RAM at the same time it sorts the
    // BAM, so give the sort the task memory minus a 40 GB genome reserve (floored
    // at 8 GB). Falls back to ~32 GB when task.memory is unset.
    def sort_ram       = task.memory ? Math.max(task.memory.toBytes() - 42949672960L, 8589934592L) : 34359738368L
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
        --limitBAMsortRAM ${sort_ram} \\
        ${quant_mode_arg}--outSAMattributes All \\
        --outSAMstrandField intronMotif \\
        --outSAMattrRGline ID:${prefix} SM:${meta.id} PL:ILLUMINA \\
        --outReadsUnmapped Fastx \\
        --outFileNamePrefix star_out/

    # STAR already emitted this coordinate-sorted (--outSAMtype BAM
    # SortedByCoordinate above), so just take it. This used to run a full
    # `samtools sort` over an already-sorted BAM — a no-op that cost ~15-25 s of
    # serial critical path in every alignment.
    mv star_out/Aligned.sortedByCoord.out.bam ${prefix}.genome_alignment.bam
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
