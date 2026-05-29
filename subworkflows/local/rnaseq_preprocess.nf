// RNA-seq preprocess: clip → [UMI extract (SE umicollapse only)] → rRNA/tRNA
// filter. Emits the genome/transcriptome-ready unaligned reads (single channel,
// SE single file or PE [R1,R2] list) plus the clip/filter logs the stats stage
// consumes. Parallel to preprocess.nf (ribo-seq).

include { CUTADAPT_CLIP_RNASEQ }  from '../../modules/local/cutadapt_clip_rnaseq.nf'
include { UMITOOLS_EXTRACT }      from '../../modules/local/umitools_extract.nf'
include { BOWTIE2_FILTER_RNASEQ } from '../../modules/local/bowtie2_filter_rnaseq.nf'

workflow RNASEQ_PREPROCESS {
    take:
    ch_reads          // [ meta(id,lane,strand,single_end), reads ]
    ch_filter_index   // value: [ index_base, [index_files] ]

    main:
    def dedup = Utils.resolve_rnaseq_dedup_method(params)

    CUTADAPT_CLIP_RNASEQ(ch_reads)

    // dedup=umicollapse extracts UMIs before filtering (SE only — validated upstream).
    ch_filter_in = (dedup == 'umicollapse') \
        ? UMITOOLS_EXTRACT(CUTADAPT_CLIP_RNASEQ.out.reads).reads \
        : CUTADAPT_CLIP_RNASEQ.out.reads

    BOWTIE2_FILTER_RNASEQ(ch_filter_in, ch_filter_index)

    // Rejoin PE R1+R2 unaligned reads into a single [meta, reads] channel.
    ch_reads_joined = BOWTIE2_FILTER_RNASEQ.out.unaligned
        .join(BOWTIE2_FILTER_RNASEQ.out.unaligned2, remainder: true)
        .map { meta, r1, r2 -> r2 ? [meta, [r1, r2]] : [meta, r1] }

    emit:
    reads_for_genome = ch_reads_joined                  // [ meta, reads ] (SE file | PE [R1,R2])
    clip_log         = CUTADAPT_CLIP_RNASEQ.out.log     // [ meta, log ]
    filter_log       = BOWTIE2_FILTER_RNASEQ.out.log    // [ meta, log ]
}
