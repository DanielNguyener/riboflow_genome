// Quality-pass filter (MAPQ + user-supplied record filter) with per-type count files.
// Ports `genome_quality_filter` (RiboFlow.groovy:576-609) and, with
// ext.presort=true, `transcriptome_sort_and_filter` (:675-701) — the STAR
// transcriptome BAM is QNAME-sorted so it must be coord-sorted after filtering.
//
// ext.mapq: `samtools view -q`, from <route>.mapping_quality_cutoff. Kept a separate
//   typed value rather than folded into ext.filter_args because on the two GENOME
//   routes the same integer also derives the stats `unique_only` mode
//   (Utils.genome_unique_only / rnaseq_genome_unique_only), which changes the shape
//   of the stats CSV. A second -q would let the filter and the stats disagree, so the
//   startup validator rejects -q inside samtools_filter_arguments.
// ext.filter_args: the rest of the user's `samtools view` record filter, verbatim from
//   <route>.samtools_filter_arguments. Set explicitly by every selector in
//   conf/modules.config so no route inherits another's; the fallback below is only a
//   null-safety net. Re-quoted token-by-token before it reaches the shell.
// ext.emit_primary_secondary (default true): when false, skip primary/secondary
// count commands (transcriptome path only needs total).
// ext.count_unique (default false): emit MAPQ-255 unique count (genome multi mode).
process SAMTOOLS_QPASS {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.qpass.bam"), path("${prefix}.qpass.bam.bai"), emit: bam
    tuple val(meta), path("${prefix}.qpass.total.count"),   emit: total_count
    tuple val(meta), path("${prefix}.qpass.primary.count"), optional: true, emit: primary_count
    tuple val(meta), path("${prefix}.qpass.secondary.count"), optional: true, emit: secondary_count
    tuple val(meta), path("${prefix}.qpass.unique.count"),  optional: true, emit: unique_count

    script:
    prefix                  = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def presort             = task.ext.presort ?: false
    def mapq                = (task.ext.mapq != null) ? task.ext.mapq : Utils.genome_mapq(params)
    def emit_primary_second = (task.ext.emit_primary_secondary != null) ? task.ext.emit_primary_secondary : true
    def count_unique        = task.ext.count_unique ?: false
    def filter_args         = Utils.shell_quote_args(
                                  (task.ext.filter_args != null) ? task.ext.filter_args.toString()
                                                                 : Utils.genome_filter_args(params))
    def sort_threads        = Math.min(task.cpus as int, 8)
    def sort_mem            = Utils.samtools_sort_mem_per_thread_mb(task)
    // PE: count fragments (first-in-pair, -f 64) so counts stay comparable to SE.
    // This -f is on the COUNTING commands only, so a user -f in filter_args cannot
    // collide with it.
    def is_pe               = meta.single_end == false
    def frag                = is_pe ? '-f 64' : ''
    // `-b -q`, not `-bq`, so the user's arguments append cleanly after the cutoff.
    def make_bam            = presort \
        ? "samtools view -h -b -q ${mapq} ${filter_args} ${bam} | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${prefix}.qpass.bam -" \
        : "samtools view -@ ${task.cpus} -b -q ${mapq} ${filter_args} ${bam} > ${prefix}.qpass.bam"
    def ps_cmd  = emit_primary_second ? """
    samtools view -@ ${task.cpus} -c -F 2304 ${frag} ${prefix}.qpass.bam > ${prefix}.qpass.primary.count
    samtools view -@ ${task.cpus} -c -f ${is_pe ? 320 : 256}  ${prefix}.qpass.bam > ${prefix}.qpass.secondary.count
    """ : ''
    def uniq_cmd = count_unique ? "samtools view -@ ${task.cpus} -c -q 255 ${frag} ${prefix}.qpass.bam > ${prefix}.qpass.unique.count" : ''
    """
    ${make_bam}
    samtools index -@ ${task.cpus} ${prefix}.qpass.bam
    samtools view -@ ${task.cpus} -c ${frag} ${prefix}.qpass.bam > ${prefix}.qpass.total.count
    ${ps_cmd}
    ${uniq_cmd}
    """

    stub:
    prefix                  = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def emit_primary_second = (task.ext.emit_primary_secondary != null) ? task.ext.emit_primary_secondary : true
    def count_unique        = task.ext.count_unique ?: false
    def ps_cmd   = emit_primary_second ? "echo 0 > ${prefix}.qpass.primary.count; echo 0 > ${prefix}.qpass.secondary.count" : ''
    def uniq_cmd = count_unique ? "echo 0 > ${prefix}.qpass.unique.count" : ''
    """
    touch ${prefix}.qpass.bam ${prefix}.qpass.bam.bai
    echo 0 > ${prefix}.qpass.total.count
    ${ps_cmd}
    ${uniq_cmd}
    """
}
