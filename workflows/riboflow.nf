// Top-level genome (ribo-seq) workflow. Builds the input channels, runs
// preprocess → genome alignment (+ optional transcriptome dedup) → stats.

include { PREPROCESS }                from '../subworkflows/local/preprocess.nf'
include { GENOME_ALIGN }              from '../subworkflows/local/genome_align.nf'
include { STAR_TRANSCRIPTOME_DEDUP }  from '../subworkflows/local/star_transcriptome_dedup.nf'
include { TRANSCRIPTOME_ALIGN }       from '../subworkflows/local/transcriptome_align.nf'
include { ALIGNMENT_STATS }           from '../subworkflows/local/alignment_stats.nf'
include { WRITE_CORRESPONDENCE }      from '../modules/local/write_correspondence.nf'

workflow RIBOFLOW {

    main:
    // ── Build per-lane input channel ──────────────────────────────────────
    // Ribo-seq lanes are single-end: one FASTQ per lane. (RiboFlow.groovy:171-177)
    def fastq_base = (params.input?.fastq_base ?: '')
    if (fastq_base && !fastq_base.endsWith('/')) fastq_base = "${fastq_base}/"

    def input_list = []
    params.input.fastq.each { sample, lanes ->
        lanes.eachWithIndex { fq, i ->
            input_list << [[id: sample, lane: (i + 1), strand: 'F'], file("${fastq_base}${fq}")]
        }
    }
    ch_reads = Channel.fromList(input_list)

    // ── Optional input existence checks (RiboFlow.groovy:200-224) ──────────
    if (params.do_check_file_existence) {
        ['SA', 'SAindex', 'Genome', 'chrNameLength.txt'].each { f ->
            assert file("${params.input.reference.genome}/${f}").exists() :
                "Missing STAR index file: ${params.input.reference.genome}/${f}"
        }
        def filter_pref = params.input.reference.filter.replaceAll('\\*', '')
        ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2'].each { s ->
            assert file("${filter_pref}.${s}").exists() :
                "Missing bowtie2 filter index file: ${filter_pref}.${s}"
        }
        input_list.each { meta, f ->
            assert f.exists() : "Missing input FASTQ: ${f}"
        }
    }

    // Sample/lane → FASTQ correspondence table.
    ch_corr = ch_reads
        .map { meta, fq -> "${meta.id}\t${meta.lane}\t${fq}" }
        .collectFile(name: 'correspondence.txt', newLine: true)
    WRITE_CORRESPONDENCE(ch_corr)

    // ── Reference channels ─────────────────────────────────────────────────
    def filter_glob = params.input.reference.filter
    def filter_base = filter_glob.split('/')[-1].replaceAll('\\*$', '').replaceAll('\\.$', '')
    ch_filter_index = Channel.value([filter_base, files(filter_glob)])
    ch_genome_index = Channel.value(file(params.input.reference.genome))

    if (params.transcriptome?.run) {
        def tx_glob = params.input.reference.transcriptome
        def tx_base = tx_glob.split('/')[-1].replaceAll('\\*$', '').replaceAll('\\.$', '')
        ch_tx_index  = Channel.value([tx_base, files(tx_glob)])
        ch_regions   = Channel.value(file(params.input.reference.regions))
        ch_lengths   = Channel.value(file(params.input.reference.transcript_lengths))
    }

    // ── Pipeline ───────────────────────────────────────────────────────────
    PREPROCESS(ch_reads, ch_filter_index)

    GENOME_ALIGN(PREPROCESS.out.reads_for_genome, ch_genome_index)

    if (Utils.do_tx_dedup(params)) {
        STAR_TRANSCRIPTOME_DEDUP(GENOME_ALIGN.out.transcriptome_bam)
    }

    if (params.transcriptome?.run) {
        TRANSCRIPTOME_ALIGN(
            PREPROCESS.out.reads_for_genome,
            ch_tx_index, ch_regions, ch_lengths
        )
    }

    ALIGNMENT_STATS(
        PREPROCESS.out.clip_log,
        PREPROCESS.out.filter_log,
        GENOME_ALIGN.out.genome_log,
        GENOME_ALIGN.out.secondary_count,
        GENOME_ALIGN.out.qpass_counts,
        GENOME_ALIGN.out.individual_dedup_counts,
        GENOME_ALIGN.out.merged_dedup_counts,
    )

    if (params.do_rnaseq && params.containsKey('rnaseq')) {
        log.warn 'RNA-seq input detected but the RNA-seq path is NOT yet migrated to DSL2 (deferred to a later stage). Ignoring rnaseq.* params.'
    }
}
