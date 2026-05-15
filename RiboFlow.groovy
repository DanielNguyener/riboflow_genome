// vim: syntax=groovy
// -*- mode: groovy;-*-

////////////////////////////////////////////////////////////////////////////////
// General Function Definitions

String get_storedir(output_type, is_rnaseq = false) {
    def base_path = params.get('output', [:]).get('intermediates', [:]).get('base', 'intermediates')
    if (is_rnaseq) {
        base_path = base_path + '/rnaseq'
    }
    new File( base_path,
              params.get('output', [:]).get('intermediates', [:]).get(output_type, output_type) )
                            .getCanonicalPath()
}

String get_rnaseq_storedir(output_type) {
    new File( params.get('output', [:]).get('intermediates', [:]).get('base', 'intermediates') + '/rnaseq',
              params.get('output', [:]).get('intermediates', [:]).get(output_type, output_type) )
                            .getCanonicalPath()
}

String get_rnaseq_publishdir(output_type) {
    new File( params.get('output', [:]).get('output', [:]).get('base', 'output') + '/rnaseq',
              params.get('output', [:]).get('output', [:]).get(output_type, output_type) )
                            .getCanonicalPath()
}

String get_publishdir(output_type, is_rnaseq = false) {
    def base_path = params.get('output', [:]).get('output', [:]).get('base', 'output')
    if (is_rnaseq) {
        base_path = base_path + '/rnaseq'
    }
    new File( base_path,
              params.get('output', [:]).get('output', [:]).get(output_type, output_type) )
                            .getCanonicalPath()
}

String get_dedup_method(String dedup_arg, String dedup_old) {
    def valid_methods = ['position', 'umicollapse', 'none']

    dedup_param = dedup_arg.toLowerCase()

    if (dedup_param != 'none') {
        if (dedup_param in valid_methods) {
            return(dedup_param)
        }
    else {
            println('Invalid deduplication method ' + dedup_param + ' . Valid methods are: ')
            println( valid_methods.join(',') )
            System.exit(1)
    }
    }
  else {
        if ( dedup_old.toLowerCase() != 'false' ) {
            return("position")
        }
    else {
            return("none")
    }
  }
}

// Per-thread memory budget for `samtools sort`. samtools allocates roughly
// this much heap per thread, so over-threaded sorts (e.g. -@ 32 × 768 MB
// default = 24 GB) easily blow past task.memory and OOM the node.
//
// We divide by the EFFECTIVE sort-thread cap (min(cpus, 8)) — not task.cpus —
// because samtools sort plateaus around 8 threads, and our call sites cap
// `-@` at 8 to match. Dividing by the real thread count keeps the per-thread
// buffer healthy when cpus=32 but only 8 threads actually run.
//
// Floor 64M = samtools' own sanity minimum. Ceiling 768M = the upstream
// default; bigger buffers don't help on the BAMs we sort.
int samtools_sort_mem_per_thread_mb(task) {
    int sort_threads = Math.min(task.cpus as int, 8)
    int est = (int) (task.memory.toMega() * 0.7 / sort_threads)
    return Math.min(768, Math.max(64, est))
}

String build_samtools_sort_cmd(sample, index, suffix, cpus) {
    return "samtools sort -@ ${cpus} -o ${sample}.${index}.${suffix}.bam"
}

String build_samtools_index_cmd(sample, index, suffix, cpus) {
    return "samtools index -@ ${cpus} ${sample}.${index}.${suffix}.bam"
}

String build_cutadapt_cmd(args, input_file, sample, index, cpus) {
    return "cutadapt --cores=${cpus} ${args} ${input_file} 2>${sample}.${index}.clipped.log | gzip -c > ${sample}.${index}.clipped.fastq.gz"
}

String build_fastqc_cmd(sample, index, cpus) {
    return """if [ ! -f ${sample}.${index}.fastq.gz ]; then
       ln -s \$input_file ${sample}.${index}.fastq.gz
    fi
    fastqc ${sample}.${index}.fastq.gz --outdir=\$PWD -t ${cpus}"""
}

////////////////////////////////////////////////////////////////////////////////

dedup_method = get_dedup_method(params.get('dedup_method', 'none').toString(),
                                params.get('deduplicate', false).toString())

do_tx_dedup = (dedup_method == 'umicollapse' || dedup_method == 'position') &&
              params.get('star', [:]).get('output_transcriptome_bam', false)

if (!params.containsKey('input')) {
    params.input = [
        reference: [
            filter: './rf_sample_data/filter/human_rtRNA*',
            genome: './rf_sample_data/genome/star_index',
        ],
        fastq: [:]
    ]
}

if (!params.containsKey('output')) {
    params.output = [
        intermediates: [
            base: 'intermediates',
            clip: 'clip',
            log: 'log',
            transcriptome_alignment: 'transcriptome_alignment',
            filter: 'filter',
            genome_alignment: 'genome_alignment',
            bam_to_bed: 'bam_to_bed',
            quality_filter: 'quality_filter',
            alignment_ribo: 'alignment_ribo',
            bigwigs: 'bigwigs'
        ],
        output: [
            base: 'output',
            log: 'log',
            fastqc: 'fastqc',
            stats: 'stats',
            bigwigs: 'bigwigs',
            alignments: 'alignments'
        ],
        individual_lane_directory: 'individual',
        merged_lane_directory: 'merged'
    ]
}

if (!params.containsKey('do_check_file_existence')) {
    params.do_check_file_existence = true
}
if (!params.containsKey('do_fastqc')) {
    params.do_fastqc = false
}
if (!params.containsKey('do_rnaseq')) {
    params.do_rnaseq = true
}
if (!params.containsKey('umi_tools_extract_arguments')) {
    params.umi_tools_extract_arguments = ''
}
if (!params.containsKey('ribo_filter_flags')) {
    params.ribo_filter_flags = 2052
}
if (!params.containsKey('do_strand_split')) {
    params.do_strand_split = false
}

rnaseq_dedup_method = params.get('rnaseq', [:]).get('dedup_method', 'position').toString()

fastq_base = params.get('input', [:]).get('fastq_base', '')
if (! fastq_base.endsWith('/') && fastq_base != '') {
    fastq_base = "${fastq_base}/"
}

Channel.from(params.get('input', [:]).fastq.collect { k, v ->
    v.collect { z -> [k, v.indexOf(z) + 1,
                                                   file("${fastq_base}${z}")] }  })
    .flatten().collate(3).into {  INPUT_SAMPLES_EXISTENCE
                                  INPUT_SAMPLES_FASTQC
                                  INPUT_SAMPLES_CLIP
                                  INPUT_SAMPLES_LOG }

INPUT_SAMPLES_LOG.flatMap { sample, index, fastq -> "${sample}\t${index}\t${fastq}" }
    .collectFile(name: 'correspondence.txt', newLine: true)
    .set { INPUT_SAMPLES_LOG_FILES }
process write_fastq_correspondence {
    executor 'local'
    publishDir get_publishdir('stats'), mode: 'move'

    input:
    file(correspondence) from INPUT_SAMPLES_LOG_FILES

    output:
    file('index_fastq_correspondence.txt')

    """
    cat ${correspondence} > index_fastq_correspondence.txt
    """
}

////////////////////////////////////////////////////////////////////////////////
// File Existence Checks

boolean file_exists(file_path) {
    this_file = file(file_path)
    assert this_file.exists()
    return true
}

boolean star_index_exists(star_dir) {
    ['SA', 'SAindex', 'Genome', 'chrNameLength.txt'].each { f ->
        file_exists("${star_dir}/${f}")
    }
    return true
}

boolean bt2_ref_exists(bt2_ref) {
    Channel.from(['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2'])
    .map { this_suffix -> file_exists("${bt2_ref}.${this_suffix}".replaceAll('\\*', '')) }
    return true
}

if (params.do_check_file_existence) {
    INPUT_SAMPLES_EXISTENCE.map { sample, index, this_file -> file_exists(this_file) }

    bt2_ref_exists(params.get('input', [:]).reference.filter)
    star_index_exists(params.get('input', [:]).reference.genome)
}

////////////////////////////////////////////////////////////////////////////////
// PROCESSES
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// RAW_FASTQC

process raw_fastqc {
    publishDir get_publishdir('fastqc') + '/raw', mode: 'copy'

    input:
    set val(sample), val(index), file(fastq) from INPUT_SAMPLES_FASTQC

    output:
    set val(sample), file("${sample}.${index}_fastqc.html"),
                       file("${sample}.${index}_fastqc.zip") into RAW_FASTQC_OUT

  when:
  params.do_fastqc

    """
    if [ ! -f ${sample}.${index}.fastq.gz ]; then
       ln -s $fastq ${sample}.${index}.fastq.gz
    fi
    fastqc ${sample}.${index}.fastq.gz --outdir=\$PWD -t ${task.cpus}
    """
}

///////////////////////////////////////////////////////////////////////////////////////
// CLIP

process clip {
    storeDir get_storedir('clip')

    input:
    set val(sample), val(index), file(fastq) from INPUT_SAMPLES_CLIP

    output:
    set val(sample), val(index), file("${sample}.${index}.clipped.fastq.gz") into CLIP_OUT
    set val(sample), val(index), file("${sample}.${index}.clipped.log") into CLIP_LOG

    """
    cutadapt --cores=${task.cpus} ${params.clip_arguments} ${fastq} 2>${sample}.${index}.clipped.log  \
     | gzip -c  > ${sample}.${index}.clipped.fastq.gz
    """
}

///////////////////////////////////////////////////////////////////////////////////////

CLIP_OUT.into { CLIP_OUT_FASTQC; CLIP_OUT_DEDUP; CLIP_OUT_FILTER }

///////////////////////////////////////////////////////////////////////////////////////
// UMI EXTRACTION

process extract_umi_via_umi_tools {
    storeDir get_storedir('umi_tools') + '/' + params.get('output', [:]).get('merged_lane_directory', 'merged')

  input:
    set val(sample), val(index), file(fastq) from CLIP_OUT_DEDUP

  output:
    set val(sample), val(index), file("${sample}.${index}.umi_extracted.fastq.gz") into UMI_EXTRACT_OUT
    set val(sample), val(index), file("${sample}.${index}.umi_extracted.log") into UMI_EXTRACT_LOG

  when:
  dedup_method == 'umicollapse'

    """
  umi_tools extract -I ${fastq} -S ${sample}.${index}.umi_extracted.fastq.gz \
     -L ${sample}.${index}.umi_extracted.log \
     ${params.umi_tools_extract_arguments ?: ''}
  """
}

///////////////////////////////////////////////////////////////////////////////////////

if (dedup_method == 'none' || dedup_method == 'position') {
    CLIP_OUT_FILTER.set { FILTER_INPUT_FASTQ }
}
else if (dedup_method == 'umicollapse') {
    UMI_EXTRACT_OUT.set { FILTER_INPUT_FASTQ }
}

///////////////////////////////////////////////////////////////////////////////////////
// CLIPPED FASTQC

process clipped_fastqc {
    publishDir get_publishdir('fastqc') + '/clipped', mode: 'copy'

    input:
    set val(sample), val(index), file(fastq)  from CLIP_OUT_FASTQC

    output:
    set val(sample), file("${sample}.${index}.clipped_fastqc.html"),
                       file("${sample}.${index}.clipped_fastqc.zip") into CLIPPED_FASTQC_OUT

    when:
    params.do_fastqc

    """
    if [ ! -f ${sample}.${index}.clipped.fastq.gz ]; then
       ln -s $fastq ${sample}.${index}.clipped.fastq.gz
    fi
    fastqc ${sample}.${index}.clipped.fastq.gz --outdir=\$PWD -t ${task.cpus}
    """
}

///////////////////////////////////////////////////////////////////////////////////////
// FILTER

FILTER_INDEX = Channel.from([[
             params.input.reference.filter
                .split('/')[-1]
                .replaceAll('\\*$', '')
                .replaceAll('\\.$', ''),
             file(params.input.reference.filter),
            ]])

process filter {
    storeDir get_storedir('filter')

    input:
    set val(sample), val(index), file(fastq) from FILTER_INPUT_FASTQ
    set val(bowtie2_index_base), file(bowtie2_index_files) from FILTER_INDEX.first()

    output:
    set val(sample), val(index), file("${sample}.${index}.filter.bam") \
        into FILTER_BAM
    set val(sample), val(index), file("${sample}.${index}.aligned.filter.fastq.gz") \
        into FILTER_ALIGNED
    set val(sample), val(index), file("${sample}.${index}.unaligned.filter.fastq.gz") \
        into FILTER_UNALIGNED
    set val(sample), val(index), file("${sample}.${index}.filter.log") \
        into FILTER_LOG

    // bowtie2 SAM piped directly into samtools sort (sort accepts SAM on stdin,
    // so the extra `samtools view -b -` step is redundant). Index for downstream
    // tools (samtools view by RG, etc.) but skip the unused idxstats sidecar.
    //
    // We decouple alignment threads from sort threads: bowtie2 doesn't scale
    // linearly past ~16 threads on the small rRNA filter index, and samtools
    // sort plateaus around 8 threads. Capping both keeps the worst-case heap
    // footprint within task.memory even when cpus=32, so two filter forks
    // can't OOM-kill each other on a node with limited free RAM.
    script:
    def aln_threads  = Math.min(task.cpus as int, 16)
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = samtools_sort_mem_per_thread_mb(task)
    """
    set -o pipefail
    bowtie2 ${params.alignment_arguments.filter} \
            -x ${bowtie2_index_base} -q ${fastq} \
            --threads ${aln_threads} \
            --al-gz ${sample}.${index}.aligned.filter.fastq.gz \
            --un-gz ${sample}.${index}.unaligned.filter.fastq.gz \
                     2> ${sample}.${index}.filter.log \
            | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${sample}.${index}.filter.bam - \
            && samtools index -@ ${sort_threads} ${sample}.${index}.filter.bam
    """
}

FILTER_UNALIGNED.set { FILTER_UNALIGNED_GENOME }

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// GENOME ALIGNMENT
// Ribo-seq libraries in this pipeline are forward-stranded: the sequenced read
// is sense to the RPF, so reads from + strand genes align to the + strand (flag 16 NOT set).
FILTER_UNALIGNED_GENOME
    .map { sample, index, fastq -> [sample, index, fastq, 'F'] }
    .set { GENOME_INPUT_CHANNEL }

GENOME_INDEX = Channel.from([[ "star_index", file(params.input.reference.genome) ]])

process genome_alignment {
    storeDir get_storedir('genome_alignment') + '/' + params.get('output', [:]).get('individual_lane_directory', 'individual')

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

input:
    set val(sample), val(index), file(fastq), val(strand_arg) from GENOME_INPUT_CHANNEL
    set val(genome_base), file(genome_files) from GENOME_INDEX.first()

output:
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.bam") \
    into GENOME_ALIGNMENT_BAM
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.bam") optional true \
    into GENOME_TRANSCRIPTOME_BAM
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.aligned.fastq.gz") \
    into GENOME_ALIGNMENT_ALIGNED
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.unaligned.fastq.gz") \
    into GENOME_ALIGNMENT_UNALIGNED
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.log") \
    into GENOME_ALIGNMENT_LOG
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.secondary.count") \
    into GENOME_ALIGNMENT_SECONDARY_COUNT
    set val(sample), val(strand_arg) into GENOME_SAMPLE_STRAND

    script:
    def emit_tx_bam = params.get('star', [:]).get('output_transcriptome_bam', false)
    def quant_mode_arg = emit_tx_bam ? '--quantMode TranscriptomeSAM \\\n        ' : ''
    def tx_bam_cmd = emit_tx_bam \
        ? "mv star_out/Aligned.toTranscriptome.out.bam ${sample}.${index}.transcriptome_alignment.bam" \
        : ''
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = samtools_sort_mem_per_thread_mb(task)
    """
    set -o pipefail
    mkdir -p star_out
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${fastq} \\
        --readFilesCommand zcat \\
        --readFilesType Fastx \\
        --genomeDir ${genome_files} \\
        --genomeLoad NoSharedMemory \\
        ${params.star.ribo_arguments} \\
        --outSAMtype BAM SortedByCoordinate \\
        ${quant_mode_arg}--outSAMattributes All \\
        --outSAMstrandField intronMotif \\
        --outSAMattrRGline ID:${sample}.${index} SM:${sample} PL:ILLUMINA \\
        --outReadsUnmapped Fastx \\
        --outFileNamePrefix star_out/

    # Defensive re-sort: STAR's multi-threaded SortedByCoordinate output can have
    # secondary-alignment records slightly out of strict coord order with
    # --outFilterMultimapNmax > 1, which breaks samtools index. Forcing an
    # explicit single-pass sort guarantees a valid coord-sorted BAM.
    # sort_threads capped at 8 — samtools sort plateaus there and bigger
    # `-@` just multiplies the per-thread memory footprint.
    samtools sort -@ ${sort_threads} -m ${sort_mem}M \\
        -o ${sample}.${index}.genome_alignment.bam \\
        star_out/Aligned.sortedByCoord.out.bam
    rm -f star_out/Aligned.sortedByCoord.out.bam
    samtools index -@ ${sort_threads} ${sample}.${index}.genome_alignment.bam

    ${tx_bam_cmd}

    samtools fastq -@ ${task.cpus} -F 4 ${sample}.${index}.genome_alignment.bam \\
        | gzip > ${sample}.${index}.genome_alignment.aligned.fastq.gz

    if [ -f star_out/Unmapped.out.mate1 ]; then
        gzip -c star_out/Unmapped.out.mate1 > ${sample}.${index}.genome_alignment.unaligned.fastq.gz
    else
        echo -n | gzip > ${sample}.${index}.genome_alignment.unaligned.fastq.gz
    fi

    cp star_out/Log.final.out ${sample}.${index}.genome_alignment.log

    # Emit secondary-alignment count next to the log. STAR reports only primary
    # alignment numbers; multi-mappers ALSO contribute secondary records to the
    # BAM (FLAG bit 256) and stats consumers need that count separately.
    samtools view -@ ${task.cpus} -c -f 256 ${sample}.${index}.genome_alignment.bam \\
        > ${sample}.${index}.genome_alignment.secondary.count
    """
}

GENOME_ALIGNMENT_ALIGNED.set { GENOME_ALIGNMENT_ALIGNED_FASTQ_FASTQC }
GENOME_ALIGNMENT_UNALIGNED.set { GENOME_ALIGNMENT_UNALIGNED_FASTQ_FASTQC }

// GENOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* GENOME ALIGNMENT FASTQC */

process genome_aligned_fastqc {
    publishDir get_publishdir('fastqc') + '/genome_aligned', mode: 'copy'

    input:
    set val(sample), val(index), file(fastq) from GENOME_ALIGNMENT_ALIGNED_FASTQ_FASTQC

    output:
    set val(sample), file("${sample}.${index}.genome.aligned_fastqc.html"),
                       file("${sample}.${index}.genome.aligned_fastqc.zip") \
                        into GENOME_ALIGNED_FASTQC_OUT

    when:
    params.do_fastqc

    """
    if [ ! -f ${sample}.${index}.genome.aligned.fastq.gz ]; then
       ln -s ${fastq} ${sample}.${index}.genome.aligned.fastq.gz
    fi

    if [ \$(stat -L -c%s ${sample}.${index}.genome.aligned.fastq.gz) -gt 20 ]; then
        fastqc ${sample}.${index}.genome.aligned.fastq.gz --outdir=\$PWD -t ${task.cpus}
    else
        echo "File is empty, skipping FastQC"
        touch ${sample}.${index}.genome.aligned_fastqc.html
        touch ${sample}.${index}.genome.aligned_fastqc.zip
    fi
    """
}

process genome_unaligned_fastqc {
    publishDir get_publishdir('fastqc') + '/genome_unaligned', mode: 'copy'

    input:
    set val(sample), val(index), file(fastq) from GENOME_ALIGNMENT_UNALIGNED_FASTQ_FASTQC

    output:
    set val(sample), file("${sample}.${index}.genome.unaligned_fastqc.html"),
                       file("${sample}.${index}.genome.unaligned_fastqc.zip") \
                        into GENOME_UNALIGNED_FASTQC_OUT

    when:
    params.do_fastqc

    """
    if [ ! -f ${sample}.${index}.genome.unaligned.fastq.gz ]; then
       ln -s ${fastq} ${sample}.${index}.genome.unaligned.fastq.gz
    fi

    if [ \$(stat -L -c%s ${sample}.${index}.genome.unaligned.fastq.gz) -gt 20 ]; then
        fastqc ${sample}.${index}.genome.unaligned.fastq.gz --outdir=\$PWD -t ${task.cpus}
    else
        echo "File is empty, skipping FastQC"
        touch ${sample}.${index}.genome.unaligned_fastqc.html
        touch ${sample}.${index}.genome.unaligned_fastqc.zip
    fi
    """
}

// GENOME ALIGNMENT FASTQC
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Genome-alignment log only feeds the stats pipeline; the merge-side and
// table-side splits had no consumers after the cleanup so they are gone.
GENOME_ALIGNMENT_LOG.set { GENOME_ALIGNMENT_LOG_STATS }

// Per-lane genome BAM only feeds genome_quality_filter (qpass). The bed
// conversion now happens AFTER qpass (see GENOME_ALIGNMENT_QPASS_BAM split
// below), and the unfiltered merge process has been removed because its
// output was unused.
GENOME_ALIGNMENT_BAM.set { GENOME_ALIGNMENT_BAM_FOR_QPASS }

if (do_tx_dedup) {
    GENOME_TRANSCRIPTOME_BAM.into {
        TRANSCRIPTOME_BAM_FOR_SORT_FILTER
        TRANSCRIPTOME_BAM_FOR_LANE_METADATA
    }
} else {
    Channel.empty().set { TRANSCRIPTOME_BAM_FOR_SORT_FILTER }
    Channel.empty().set { TRANSCRIPTOME_BAM_FOR_LANE_METADATA }
}

process genome_quality_filter {
    storeDir get_storedir('quality_filter')

    input:
        set val(sample), val(index), file(bam) from GENOME_ALIGNMENT_BAM_FOR_QPASS

    output:
        set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.bam") \
            into GENOME_ALIGNMENT_QPASS_BAM
        set val(sample), val(index),
            file("${sample}.${index}.qpass.total.count"),
            file("${sample}.${index}.qpass.primary.count"),
            file("${sample}.${index}.qpass.secondary.count") \
            into GENOME_QPASS_COUNTS_OUT

    // -F ribo_filter_flags: drop unmapped (4) + supplementary (2048); KEEP secondary (256) so
    // multi-mapper alignments contribute to every locus in the final bigwig.
    // Three counts feed the per-step alignment breakdown in the stats CSV:
    //   total      = all records in qpass BAM
    //   primary    = -F 2304 (drops 256 secondary + 2048 supplementary; one record per read)
    //   secondary  = -f 256
    """
    samtools view -@ ${task.cpus} -bq ${params.mapping_quality_cutoff} \\
        -F ${params.ribo_filter_flags} ${bam} \\
        > ${sample}.${index}.genome_alignment.qpass.bam
    samtools index -@ ${task.cpus} ${sample}.${index}.genome_alignment.qpass.bam
    samtools view -@ ${task.cpus} -c        ${sample}.${index}.genome_alignment.qpass.bam \\
        > ${sample}.${index}.qpass.total.count
    samtools view -@ ${task.cpus} -c -F 2304 ${sample}.${index}.genome_alignment.qpass.bam \\
        > ${sample}.${index}.qpass.primary.count
    samtools view -@ ${task.cpus} -c -f 256  ${sample}.${index}.genome_alignment.qpass.bam \\
        > ${sample}.${index}.qpass.secondary.count
    """
}

// QPASS count tuples carry (total, primary, secondary) for the new stats schema.
GENOME_QPASS_COUNTS_OUT.into {
    GENOME_QPASS_COUNTS
    GENOME_QPASS_COUNTS_FOR_STATS
    GENOME_QPASS_COUNTS_FOR_NONE_DEDUP
}

GENOME_ALIGNMENT_QPASS_BAM.into {
    GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE
    GENOME_ALIGNMENT_QPASS_BAM_FOR_BED        // per-lane BED conversion (replaces
                                              // raw-BAM conversion; matches the
                                              // transcriptome path)
    GENOME_ALIGNMENT_QPASS_BAM_FOR_NODEDUP_PUB // per-lane qpass BAM publish on the
                                              // dedup_method == 'none' path
}
GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE.map { sample, index, bam -> [sample, bam] }.groupTuple()
.set { GENOME_ALIGNMENT_QPASS_GROUPED_BAM }

// Single clean merge process — output is the merged qpass BAM regardless of
// dedup method. When dedup_method == 'none' the BAM is also published as a
// user-facing artifact (output/alignments/ribo/merged/) because it IS the
// final BAM. When dedup is on, the BAM stays in intermediates and feeds the
// dedup stages.
process merge_genome_qpass_bam {
    storeDir get_storedir('genome_alignment') + '/' + params.output.merged_lane_directory
    publishDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory,
               mode: 'copy', enabled: dedup_method == 'none'

    input:
        set val(sample), file(bam_files) from GENOME_ALIGNMENT_QPASS_GROUPED_BAM

    output:
        set val(sample), file("${sample}.genome.qpass.merged.bam"),
                         file("${sample}.genome.qpass.merged.bam.bai") \
            into GENOME_MERGED_QPASS_BAM

    """
    samtools merge -@ ${task.cpus} ${sample}.genome.qpass.merged.bam ${bam_files}
    samtools index -@ ${task.cpus} ${sample}.genome.qpass.merged.bam
    """
}

// Route the merged qpass BAM by dedup method. For 'none' it feeds the bigwig
// stage AND a new merged-qpass-bed publish step; for 'umicollapse'/'position'
// it feeds the dedup chain.
if (dedup_method == 'none') {
    GENOME_MERGED_QPASS_BAM.set { GENOME_MERGED_BAM_QPASS_NONE }
    Channel.empty().set { GENOME_MERGED_BAM_FOR_UMICOLLAPSE_DEDUP }
    Channel.empty().set { GENOME_MERGED_BAM_FOR_POSITION_BAM_CONV }
} else {
    GENOME_MERGED_QPASS_BAM.into {
        GENOME_MERGED_BAM_FOR_UMICOLLAPSE_DEDUP
        GENOME_MERGED_BAM_FOR_POSITION_BAM_CONV
    }
    Channel.empty().set { GENOME_MERGED_BAM_QPASS_NONE }
}

///////////////////////////////////////////////////////////////////////////////
//   T R A N S C R I P T O M E   D E D U P L I C A T I O N
//   Processes transcriptome-projected BAMs (star.output_transcriptome_bam: true)
//   through the same dedup pipeline as genome BAMs.
///////////////////////////////////////////////////////////////////////////////

// Coordinate-sort the QNAME-sorted STAR transcriptome BAM and apply qpass filter.
process transcriptome_sort_and_filter {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.individual_lane_directory

    input:
    set val(sample), val(index), file(bam) from TRANSCRIPTOME_BAM_FOR_SORT_FILTER

    output:
    set val(sample), val(index),
        file("${sample}.${index}.transcriptome_alignment.qpass.bam"),
        file("${sample}.${index}.transcriptome_alignment.qpass.bam.bai") \
        into TRANSCRIPTOME_QPASS_BAM

    when:
    do_tx_dedup

    // STAR emits transcriptome BAM in QNAME order. Quality-filter on SAM then
    // coord-sort, all in one pipe — no intermediate `sorted.bam` on disk.
    // sort_threads capped at 8 (samtools sort scales poorly past that).
    script:
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = samtools_sort_mem_per_thread_mb(task)
    """
    samtools view -h -bq ${params.mapping_quality_cutoff} -F ${params.ribo_filter_flags} ${bam} \\
        | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${sample}.${index}.transcriptome_alignment.qpass.bam -
    samtools index -@ ${sort_threads} ${sample}.${index}.transcriptome_alignment.qpass.bam
    """
}

if (do_tx_dedup) {
    TRANSCRIPTOME_QPASS_BAM.into {
        TRANSCRIPTOME_QPASS_BAM_FOR_BED    // position path: per-lane BED conversion
        TRANSCRIPTOME_QPASS_BAM_FOR_MERGE  // both paths: lane merging
    }
} else {
    Channel.empty().set { TRANSCRIPTOME_QPASS_BAM_FOR_BED }
    Channel.empty().set { TRANSCRIPTOME_QPASS_BAM_FOR_MERGE }
}

// ── POSITION DEDUP PATH ───────────────────────────────────────────────────────

process transcriptome_qpass_bam_to_bed {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.individual_lane_directory

    input:
    set val(sample), val(index), file(bam), file(bai) from TRANSCRIPTOME_QPASS_BAM_FOR_BED

    output:
    set val(sample), val(index),
        file("${sample}.${index}.transcriptome.with_sample_index.bed") \
        into TRANSCRIPTOME_BED_WITH_INDEX

    when:
    dedup_method == 'position' && do_tx_dedup

    """
    bamToBed -i ${bam} \
        | awk -v s=${sample}.${index} '{ print \$0"\t"s }' \
        > ${sample}.${index}.transcriptome.with_sample_index.bed
    """
}

if (dedup_method == 'position' && do_tx_dedup) {
    TRANSCRIPTOME_BED_WITH_INDEX
        .map { sample, index, bed -> [sample, bed] }
        .groupTuple()
        .set { TRANSCRIPTOME_BED_WITH_INDEX_GROUPED }
} else {
    Channel.empty().set { TRANSCRIPTOME_BED_WITH_INDEX_GROUPED }
}

process merge_transcriptome_bed_for_position_dedup {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bed_files) from TRANSCRIPTOME_BED_WITH_INDEX_GROUPED

    output:
    set val(sample), file("${sample}.transcriptome.merged.pre_dedup.bed") \
        into TRANSCRIPTOME_MERGED_BED_PRE_DEDUP

    when:
    dedup_method == 'position' && do_tx_dedup

    """
    cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n \
        > ${sample}.transcriptome.merged.pre_dedup.bed
    """
}

process transcriptome_deduplicate_position {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.merged_lane_directory

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

    input:
    set val(sample), file(bed) from TRANSCRIPTOME_MERGED_BED_PRE_DEDUP

    output:
    set val(sample), file("${sample}.transcriptome.post_dedup.bed") \
        into TRANSCRIPTOME_POST_DEDUP_BED_RAW

    when:
    dedup_method == 'position' && do_tx_dedup

    """
    sort -k1,1 -k2,2n -k3,3n ${bed} > sorted.bed
    rfc dedup -i sorted.bed -o ${sample}.transcriptome.post_dedup.bed
    rm sorted.bed
    """
}

if (dedup_method == 'position' && do_tx_dedup) {
    TRANSCRIPTOME_POST_DEDUP_BED_RAW.into {
        TRANSCRIPTOME_POST_DEDUP_BED_FOR_SEPARATION
        TRANSCRIPTOME_POST_DEDUP_BED_FOR_BAM_CONV
    }

    TRANSCRIPTOME_BAM_FOR_LANE_METADATA
        .map { sample, index, bam -> [sample, index] }
        .combine(TRANSCRIPTOME_POST_DEDUP_BED_FOR_SEPARATION, by: 0)
        .set { TRANSCRIPTOME_BED_FOR_SEPARATION }
} else {
    Channel.empty().set { TRANSCRIPTOME_POST_DEDUP_BED_FOR_BAM_CONV }
    Channel.empty().set { TRANSCRIPTOME_BED_FOR_SEPARATION }
}

process separate_transcriptome_bed_post_dedup {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.individual_lane_directory

    input:
    set val(sample), val(index), file(bed) from TRANSCRIPTOME_BED_FOR_SEPARATION

    output:
    set val(sample), val(index),
        file("${sample}.${index}.transcriptome.post_dedup.bed") \
        into TRANSCRIPTOME_INDIVIDUAL_POST_DEDUP_BED

    when:
    dedup_method == 'position' && do_tx_dedup

    """
    awk -v s=${sample}.${index} '\$7 == s { print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6 }' \
        ${bed} > ${sample}.${index}.transcriptome.post_dedup.bed
    """
}

if (dedup_method == 'position' && do_tx_dedup) {
    TRANSCRIPTOME_QPASS_BAM_FOR_MERGE
        .map { sample, index, bam, bai -> [sample, bam] }
        .groupTuple()
        .set { TRANSCRIPTOME_QPASS_GROUPED_FOR_POSITION }
} else {
    Channel.empty().set { TRANSCRIPTOME_QPASS_GROUPED_FOR_POSITION }
}

process merge_transcriptome_qpass_bam_for_position {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam_files) from TRANSCRIPTOME_QPASS_GROUPED_FOR_POSITION

    output:
    set val(sample), file("${sample}.transcriptome.qpass.merged.bam"),
                     file("${sample}.transcriptome.qpass.merged.bam.bai") \
        into TRANSCRIPTOME_MERGED_QPASS_BAM_FOR_CONV

    when:
    dedup_method == 'position' && do_tx_dedup

    """
    samtools merge -@ ${task.cpus} ${sample}.transcriptome.qpass.merged.bam ${bam_files}
    samtools index -@ ${task.cpus} ${sample}.transcriptome.qpass.merged.bam
    """
}

if (dedup_method == 'position' && do_tx_dedup) {
    TRANSCRIPTOME_POST_DEDUP_BED_FOR_BAM_CONV
        .join(TRANSCRIPTOME_MERGED_QPASS_BAM_FOR_CONV)
        .map { sample, bed, bam, bai -> [sample, bed, bam] }
        .set { TRANSCRIPTOME_BED_BAM_FOR_CONV }
} else {
    Channel.empty().set { TRANSCRIPTOME_BED_BAM_FOR_CONV }
}

process transcriptome_convert_dedup_bed_to_bam_position {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(dedup_bed), file(qpass_bam) from TRANSCRIPTOME_BED_BAM_FOR_CONV

    output:
    set val(sample), file("${sample}.transcriptome.post_dedup.bam"),
                     file("${sample}.transcriptome.post_dedup.bam.bai") \
        into TRANSCRIPTOME_MERGED_DEDUP_BAM_POSITION

    when:
    dedup_method == 'position' && do_tx_dedup

    """
    set +u
    source \$(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome
    set -u
    rfc extract-dedup-reads \\
        --bam ${qpass_bam} \\
        --bed ${dedup_bed} \\
        --output ${sample}.transcriptome.post_dedup.bam
    samtools index -@ ${task.cpus} ${sample}.transcriptome.post_dedup.bam
    """
}

// ── UMICOLLAPSE DEDUP PATH ────────────────────────────────────────────────────

if (dedup_method == 'umicollapse' && do_tx_dedup) {
    TRANSCRIPTOME_QPASS_BAM_FOR_MERGE
        .map { sample, index, bam, bai -> [sample, bam] }
        .groupTuple()
        .set { TRANSCRIPTOME_QPASS_GROUPED_FOR_UMI }
} else {
    Channel.empty().set { TRANSCRIPTOME_QPASS_GROUPED_FOR_UMI }
}

process merge_transcriptome_bam_for_umicollapse {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam_files) from TRANSCRIPTOME_QPASS_GROUPED_FOR_UMI

    output:
    set val(sample), file("${sample}.transcriptome.qpass.merged.bam"),
                     file("${sample}.transcriptome.qpass.merged.bam.bai") \
        into TRANSCRIPTOME_MERGED_BAM_FOR_UMI_DEDUP

    when:
    dedup_method == 'umicollapse' && do_tx_dedup

    """
    samtools merge -@ ${task.cpus} ${sample}.transcriptome.qpass.merged.bam ${bam_files}
    samtools index -@ ${task.cpus} ${sample}.transcriptome.qpass.merged.bam
    """
}

process transcriptome_deduplicate_umicollapse {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam), file(bai) from TRANSCRIPTOME_MERGED_BAM_FOR_UMI_DEDUP

    output:
    set val(sample), file("${sample}.transcriptome.dedup.bam") \
        into TRANSCRIPTOME_UMI_DEDUP_BAM

    when:
    dedup_method == 'umicollapse' && do_tx_dedup

    """
    java11 -Xms512m -Xmx32g -Xss256m \\
        umicollapse.main.Main bam \\
        -i ${bam} \\
        -o ${sample}.transcriptome.dedup.bam \\
        --umi-sep "_" \\
        --algo dir \\
        --merge mapqual \\
        --two-pass \\
        ${params.get('umicollapse_arguments', '')}
    samtools index -@ ${task.cpus} ${sample}.transcriptome.dedup.bam
    """
}

if (dedup_method == 'umicollapse' && do_tx_dedup) {
    TRANSCRIPTOME_UMI_DEDUP_BAM
        .cross(TRANSCRIPTOME_BAM_FOR_LANE_METADATA.map { sample, index, bam -> [sample, index] })
        .map { merged, individual -> [individual[0], individual[1], merged[1]] }
        .set { TRANSCRIPTOME_UMI_DEDUP_BAM_FOR_SPLITTING }
} else {
    Channel.empty().set { TRANSCRIPTOME_UMI_DEDUP_BAM_FOR_SPLITTING }
}

process split_transcriptome_dedup_bam_to_individual {
    storeDir get_storedir('transcriptome_alignment') + '/' + params.output.individual_lane_directory

    input:
    set val(sample), val(index), file(merged_bam) from TRANSCRIPTOME_UMI_DEDUP_BAM_FOR_SPLITTING

    output:
    set val(sample), val(index),
        file("${sample}.${index}.transcriptome.post_dedup.bam"),
        file("${sample}.${index}.transcriptome.post_dedup.bam.bai") \
        into TRANSCRIPTOME_INDIVIDUAL_DEDUP_BAM_UMI

    when:
    dedup_method == 'umicollapse' && do_tx_dedup

    """
    if ! samtools view -H ${merged_bam} | grep -q "^@RG"; then
        echo "ERROR: No read groups found in ${merged_bam}" >&2; exit 1
    fi
    samtools view -@ ${task.cpus} -B -r ${sample}.${index} ${merged_bam} \\
        -o ${sample}.${index}.transcriptome.post_dedup.bam
    samtools index -@ ${task.cpus} ${sample}.${index}.transcriptome.post_dedup.bam
    """
}

///////////////////////////////////////////////////////////////////////////////

// Per-lane qpass BAM → BED. Matches the transcriptome path (bed derived from
// the qpass-filtered BAM, not the raw STAR output). When dedup_method == 'none'
// these BEDs are the final per-lane artifacts and get published; otherwise
// they only feed the position-dedup chain (and a sample/lane lookup helper).
process individual_genome_bam_to_bed {
    storeDir get_storedir('bam_to_bed') + '/' + params.output.individual_lane_directory
    publishDir get_publishdir('alignments') + '/ribo/' + params.output.individual_lane_directory,
               mode: 'copy', enabled: dedup_method == 'none'

    input:
        set val(sample), val(index), file(bam) from GENOME_ALIGNMENT_QPASS_BAM_FOR_BED

    output:
        set val(sample), val(index), file("${sample}.${index}.genome.qpass.bed") \
            into GENOME_INDIVIDUAL_QPASS_BED

    """
    if [ \$(samtools view -@ ${task.cpus} -c ${bam}) -eq 0 ]; then
        touch ${sample}.${index}.genome.qpass.bed
    else
        bamToBed -i ${bam} > ${sample}.${index}.genome.qpass.bed
    fi
    """
}

GENOME_INDIVIDUAL_QPASS_BED.into {
    GENOME_INDIVIDUAL_QPASS_BED_FOR_INDEX     // position-dedup: add sample-index column
    GENOME_INDIVIDUAL_QPASS_BED_FOR_LOOKUP    // sample/index lookup helper (umicollapse split routing)
    GENOME_INDIVIDUAL_QPASS_BED_FOR_NONE_MERGE // dedup_method='none' → concat for merged qpass BED publish
}

// Position-dedup-only: tag each BED record with the sample/lane identifier so
// the rfc-dedup output can be split back per-lane via awk on column 7.
process add_sample_index_col_to_genome_bed {
    storeDir get_storedir('bam_to_bed') + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bed) from GENOME_INDIVIDUAL_QPASS_BED_FOR_INDEX

    output:
        set val(sample), file("${sample}.${index}.genome.with_sample_index.bed") \
         into GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED

    when:
        dedup_method == 'position'

        """
    awk -v this_sample=${sample}.${index} \
     '{ print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"this_sample) }' ${bed} \
          > ${sample}.${index}.genome.with_sample_index.bed
    """
}
///////////////////////////////////////////////////////////////////////////////
//          G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////
//
// The previous `merge_genome_alignment` (per-lane unfiltered BAM merge +
// FASTQ/log/CSV concat) and `genome_bam_to_bed` (BED from that merged BAM)
// have been removed: their outputs were never consumed. Position-dedup now
// builds its pre-dedup BED from the per-lane qpass BAM via
// `individual_genome_bam_to_bed` → `add_sample_index_col_to_genome_bed` →
// `merge_genome_bed_for_position_dedup`, matching the transcriptome path.
//
if (dedup_method == 'position') {
    GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED
        .map { sample, bed -> [sample, bed] }
        .groupTuple()
        .set { GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED }

    process merge_genome_bed_for_position_dedup {
        storeDir get_storedir('bam_to_bed') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(bed_files) from GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED

        output:
            set val(sample), file("${sample}.genome.merged.pre_dedup.bed") \
                into GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION

        when:
        dedup_method == 'position'

        """
        cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.genome.merged.pre_dedup.bed
        """
    }
} else {
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION }
}

process genome_deduplicate_position {
    storeDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

    input:
        set val(sample), file(bed) from GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION

    output:
        set val(sample), file("${sample}.genome.post_dedup.bed") \
         into GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW

    when:
    dedup_method == 'position'

        """
    sort -k1,1 -k2,2n -k3,3n ${bed} > ${sample}.genome.sorted.bed
    rfc dedup -i ${sample}.genome.sorted.bed -o ${sample}.genome.post_dedup.bed
    """
}

///////////////////////////////////////////////////////////////////////////////
// Dedup post-processing — produces (a) per-lane and merged post-dedup BED/BAM
// final artifacts, and (b) per-lane + sample-level alignment count tuples
// (total/primary/secondary) consumed by the stats pipeline.
//
// The previous module had ~7 small bookkeeping processes that re-derived BED
// counts via `wc -l` and re-derived BEDs via a second `bamToBed`. Both are
// gone: counts now come from `samtools view -c` on the actual BAM, and the
// merged umicollapse BED is concatenated from the per-lane BEDs that the
// split step already emits.
///////////////////////////////////////////////////////////////////////////////

if (dedup_method == 'position') {

    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.into {
        GENOME_BED_FOR_DEDUP_POSITION_FOR_SEPARATION
        GENOME_BED_FOR_DEDUP_POSITION_FOR_BAM
    }

    // Awk-split the merged post-dedup BED on the sample-index column. Per-lane
    // count tuple is emitted from the per-lane BED:
    //   total   = wc -l (one record per alignment, primary or secondary)
    //   primary = unique read names in col 4 (each read has exactly one primary)
    //   secondary = total - primary
    GENOME_QPASS_COUNTS_FOR_STATS
        .map { sample, index, total, primary, secondary -> [sample, index] }
        .unique()
        .combine(GENOME_BED_FOR_DEDUP_POSITION_FOR_SEPARATION, by: 0)
        .set { GENOME_BED_FOR_POSITION_SEPARATION }

    process separate_genome_bed_post_dedup {
        storeDir get_publishdir('alignments') + '/ribo/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(bed) from GENOME_BED_FOR_POSITION_SEPARATION

        output:
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") \
                into GENOME_INDIVIDUAL_POST_DEDUP_BED_POSITION
            set val(sample), val(index),
                file("${sample}.${index}.dedup.total.count"),
                file("${sample}.${index}.dedup.primary.count"),
                file("${sample}.${index}.dedup.secondary.count") \
                into GENOME_INDIVIDUAL_DEDUP_COUNT_POSITION

        when:
        dedup_method == 'position'

        """
        awk -v s=${sample}.${index} \\
            '\$7 == s { print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6 }' \\
            ${bed} > ${sample}.${index}.genome.post_dedup.bed
        total=\$(wc -l < ${sample}.${index}.genome.post_dedup.bed)
        primary=\$(awk '{print \$4}' ${sample}.${index}.genome.post_dedup.bed | sort -u | wc -l)
        secondary=\$((total - primary))
        echo \${total}     > ${sample}.${index}.dedup.total.count
        echo \${primary}   > ${sample}.${index}.dedup.primary.count
        echo \${secondary} > ${sample}.${index}.dedup.secondary.count
        """
    }

    // Extract reads from the qpass merged BAM that match the rfc-dedup BED →
    // sample-level post-dedup BAM (the bigwig source). Counts come from the
    // resulting BAM, not from the BED, so secondary alignments are tracked
    // correctly even though rfc dedup operates on coordinate-only BED.
    GENOME_BED_FOR_DEDUP_POSITION_FOR_BAM
        .join(GENOME_MERGED_BAM_FOR_POSITION_BAM_CONV)
        .map { sample, dedup_bed, bam, bai -> [sample, dedup_bed, bam] }
        .set { GENOME_BED_FOR_BAM_CONVERSION_POSITION }

    process genome_convert_dedup_bed_to_bam_position {
        storeDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory

        input:
        set val(sample), file(dedup_bed), file(qpass_bam) from GENOME_BED_FOR_BAM_CONVERSION_POSITION

        output:
        set val(sample), file("${sample}.post_dedup.bam"), file("${sample}.post_dedup.bam.bai") \
            into GENOME_POST_DEDUP_BAM_POSITION_MERGED
        set val(sample),
            file("${sample}.merged_dedup.total.count"),
            file("${sample}.merged_dedup.primary.count"),
            file("${sample}.merged_dedup.secondary.count") \
            into GENOME_MERGED_DEDUP_COUNTS_POSITION

        when:
        dedup_method == 'position'

        """
        set +u
        source \$(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome
        set -u
        rfc extract-dedup-reads \\
            --bam ${qpass_bam} \\
            --bed ${dedup_bed} \\
            --output ${sample}.post_dedup.bam
        samtools index -@ ${task.cpus} ${sample}.post_dedup.bam
        samtools view -@ ${task.cpus} -c        ${sample}.post_dedup.bam \\
            > ${sample}.merged_dedup.total.count
        samtools view -@ ${task.cpus} -c -F 2304 ${sample}.post_dedup.bam \\
            > ${sample}.merged_dedup.primary.count
        samtools view -@ ${task.cpus} -c -f 256  ${sample}.post_dedup.bam \\
            > ${sample}.merged_dedup.secondary.count
        """
    }

    GENOME_POST_DEDUP_BAM_POSITION_MERGED.set { GENOME_BAM_FOR_BIGWIG_FINAL }
    GENOME_INDIVIDUAL_DEDUP_COUNT_POSITION.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    GENOME_MERGED_DEDUP_COUNTS_POSITION.set { GENOME_MERGED_DEDUP_COUNTS }

} else if (dedup_method == 'umicollapse') {

    // umicollapse: BAM-in, BAM-out. Emit BAM+BAI as one tuple (no follow-up
    // index step needed) and the three alignment counts in the same process.
    process genome_deduplicate_umicollapse {
        storeDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory

        input:
        set val(sample), file(bam), file(bai) from GENOME_MERGED_BAM_FOR_UMICOLLAPSE_DEDUP

        output:
        set val(sample), file("${sample}.dedup.bam"), file("${sample}.dedup.bam.bai") \
            into GENOME_UMI_DEDUP_BAM
        set val(sample),
            file("${sample}.merged_dedup.total.count"),
            file("${sample}.merged_dedup.primary.count"),
            file("${sample}.merged_dedup.secondary.count") \
            into GENOME_MERGED_DEDUP_COUNTS_UMI

        when:
        dedup_method == 'umicollapse'

        """
        java11 -Xms512m -Xmx32g -Xss256m \\
            umicollapse.main.Main bam \\
            -i ${bam} \\
            -o ${sample}.dedup.bam \\
            --umi-sep "_" \\
            --algo dir \\
            --merge mapqual \\
            --two-pass \\
            ${params.get('umicollapse_arguments', '')}
        samtools index -@ ${task.cpus} ${sample}.dedup.bam
        samtools view -@ ${task.cpus} -c        ${sample}.dedup.bam \\
            > ${sample}.merged_dedup.total.count
        samtools view -@ ${task.cpus} -c -F 2304 ${sample}.dedup.bam \\
            > ${sample}.merged_dedup.primary.count
        samtools view -@ ${task.cpus} -c -f 256  ${sample}.dedup.bam \\
            > ${sample}.merged_dedup.secondary.count
        """
    }

    GENOME_UMI_DEDUP_BAM.into {
        GENOME_UMI_DEDUP_BAM_FOR_BIGWIG
        GENOME_UMI_DEDUP_BAM_FOR_SPLITTING
    }
    GENOME_UMI_DEDUP_BAM_FOR_BIGWIG.set { GENOME_BAM_FOR_BIGWIG_FINAL }

    // Per-lane split: BAM is split by @RG (needs the merged BAM only; bai
    // sidecar from the dedup step is already in scope). bamToBed runs once
    // per lane — no separate bamToBed pass on the merged BAM is needed.
    GENOME_UMI_DEDUP_BAM_FOR_SPLITTING
        .cross(GENOME_INDIVIDUAL_QPASS_BED_FOR_LOOKUP.map { sample, index, bed -> [sample, index] }.unique())
        .map { merged_data, individual_data ->
            [individual_data[0], individual_data[1], merged_data[1], merged_data[2]]
        }
        .set { GENOME_DEDUP_BAM_FOR_SPLITTING }

    process split_genome_dedup_bam_to_individual {
        storeDir get_publishdir('alignments') + '/ribo/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(merged_bam), file(merged_bai) from GENOME_DEDUP_BAM_FOR_SPLITTING

        output:
            set val(sample), val(index),
                file("${sample}.${index}.post_dedup.bam"),
                file("${sample}.${index}.post_dedup.bam.bai") \
                into GENOME_INDIVIDUAL_DEDUP_BAM_UMI
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") \
                into GENOME_INDIVIDUAL_POST_DEDUP_BED_UMI
            set val(sample), val(index),
                file("${sample}.${index}.dedup.total.count"),
                file("${sample}.${index}.dedup.primary.count"),
                file("${sample}.${index}.dedup.secondary.count") \
                into GENOME_INDIVIDUAL_DEDUP_COUNT_UMI

        when:
        dedup_method == 'umicollapse'

        """
        if ! samtools view -H ${merged_bam} | grep -q "^@RG"; then
            echo "ERROR: No read groups found in merged BAM ${merged_bam}" >&2
            exit 1
        fi
        samtools view -@ ${task.cpus} -B -r ${sample}.${index} ${merged_bam} \\
            -o ${sample}.${index}.post_dedup.bam
        samtools index -@ ${task.cpus} ${sample}.${index}.post_dedup.bam
        samtools view -@ ${task.cpus} -c        ${sample}.${index}.post_dedup.bam \\
            > ${sample}.${index}.dedup.total.count
        samtools view -@ ${task.cpus} -c -F 2304 ${sample}.${index}.post_dedup.bam \\
            > ${sample}.${index}.dedup.primary.count
        samtools view -@ ${task.cpus} -c -f 256  ${sample}.${index}.post_dedup.bam \\
            > ${sample}.${index}.dedup.secondary.count
        if [ \$(cat ${sample}.${index}.dedup.total.count) -eq 0 ]; then
            touch ${sample}.${index}.genome.post_dedup.bed
        else
            bamToBed -i ${sample}.${index}.post_dedup.bam \\
                > ${sample}.${index}.genome.post_dedup.bed
        fi
        """
    }

    // Concat per-lane post-dedup BEDs → merged post-dedup BED (publish-only;
    // no downstream channel consumer, replaces the deleted bamToBed-on-merged
    // pass that produced the same artifact at higher cost).
    GENOME_INDIVIDUAL_POST_DEDUP_BED_UMI
        .map { sample, index, bed -> [sample, bed] }
        .groupTuple()
        .set { GENOME_INDIVIDUAL_POST_DEDUP_BED_UMI_GROUPED }

    process concat_genome_post_dedup_bed_umi {
        storeDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory

        input:
            set val(sample), file(bed_files) from GENOME_INDIVIDUAL_POST_DEDUP_BED_UMI_GROUPED

        output:
            set val(sample), file("${sample}.genome.post_dedup.bed") \
                into GENOME_MERGED_POST_DEDUP_BED_UMI

        when:
        dedup_method == 'umicollapse'

        """
        cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.genome.post_dedup.bed
        """
    }

    GENOME_INDIVIDUAL_DEDUP_COUNT_UMI.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    GENOME_MERGED_DEDUP_COUNTS_UMI.set { GENOME_MERGED_DEDUP_COUNTS }

} else {
    // dedup_method == 'none' — qpass IS the final output. Bigwig comes from
    // the merged qpass BAM. Per-lane and merged dedup counts mirror qpass.
    GENOME_MERGED_BAM_QPASS_NONE.set { GENOME_BAM_FOR_BIGWIG_FINAL }
    GENOME_QPASS_COUNTS_FOR_NONE_DEDUP.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    // No sample-level merged dedup count override needed — the merged-stats
    // helper falls back to summing per-lane qpass counts when this is empty.
    Channel.empty().set { GENOME_MERGED_DEDUP_COUNTS }

    // For 'none' we also publish the merged qpass BED (concat of per-lane
    // qpass BEDs) under output/alignments/ribo/merged/.
    GENOME_INDIVIDUAL_QPASS_BED_FOR_NONE_MERGE
        .map { sample, index, bed -> [sample, bed] }
        .groupTuple()
        .set { GENOME_INDIVIDUAL_QPASS_BED_GROUPED_FOR_NONE }

    process concat_genome_qpass_bed_for_publish_none {
        storeDir get_storedir('bam_to_bed') + '/' + params.output.merged_lane_directory
        publishDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory,
                   mode: 'copy'

        input:
            set val(sample), file(bed_files) from GENOME_INDIVIDUAL_QPASS_BED_GROUPED_FOR_NONE

        output:
            set val(sample), file("${sample}.genome.qpass.bed") \
                into GENOME_MERGED_QPASS_BED_NONE

        when:
        dedup_method == 'none'

        """
        cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.genome.qpass.bed
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
//  R I B O - S E Q   B I G W I G   G E N E R A T I O N
///////////////////////////////////////////////////////////////////////////////

GENOME_SAMPLE_STRAND
    .groupTuple()
    .map { sample, strands -> [sample, strands[0]] }
    .set { GENOME_SAMPLE_STRAND_UNIQUE }

GENOME_BAM_FOR_BIGWIG_FINAL
    .join(GENOME_SAMPLE_STRAND_UNIQUE)
    .into { GENOME_BAM_FOR_BIGWIG_WITH_STRAND
            GENOME_BAM_FOR_STRAND_SPLIT }

process genome_create_strand_specific_bigwigs {
    storeDir get_publishdir('bigwigs') + '/ribo'

   

    input:
    set val(sample), file(bam), file(bai), val(strand_arg) from GENOME_BAM_FOR_BIGWIG_WITH_STRAND

    output:
    set val(sample), file("${sample}.ribo.plus.bigWig"), \
                     file("${sample}.ribo.minus.bigWig") \
        into GENOME_STRAND_SPECIFIC_BIGWIGS

    when:
    dedup_method == 'umicollapse' || dedup_method == 'none' || dedup_method == 'position'

    // deepTools --filterRNAstrand assumes dUTP/reverse-stranded libraries, so
    // its "forward"/"reverse" labels are inverted relative to forward-stranded
    // ribo-seq (read is sense to the RPF). Map plus<->reverse and minus<->forward.
    //
    // bw_threads capped at 8: each `bamCoverage -p N` forks N Python worker
    // processes via multiprocessing, and each fork copies the BAM index +
    // working buffers. 32 forks per BAM × 2 BAMs per sample × multiple samples
    // exhausts the kernel fork table and triggers ENOMEM on fork().
    script:
    def bw_threads = Math.min(task.cpus as int, 8)
    """
    if [ "${strand_arg}" == "F" ] || [ "${strand_arg}" == "FR" ]; then
        bamCoverage -b ${bam} -o ${sample}.ribo.plus.bigWig \
            --filterRNAstrand reverse --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
        bamCoverage -b ${bam} -o ${sample}.ribo.minus.bigWig \
            --filterRNAstrand forward --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
    else
        bamCoverage -b ${bam} -o ${sample}.ribo.plus.bigWig \
            --filterRNAstrand forward --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
        bamCoverage -b ${bam} -o ${sample}.ribo.minus.bigWig \
            --filterRNAstrand reverse --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig
    fi
    """
}

process genome_split_stranded_bam {
    storeDir get_publishdir('alignments') + '/ribo/stranded'

    input:
    set val(sample), file(bam), file(bai), val(strand_arg) from GENOME_BAM_FOR_STRAND_SPLIT

    output:
    set val(sample), file("${sample}.ribo.plus.bam"),  file("${sample}.ribo.plus.bam.bai"),
                     file("${sample}.ribo.minus.bam"), file("${sample}.ribo.minus.bam.bai") \
        into GENOME_STRANDED_SPLIT_BAMS

    when:
    params.get('do_strand_split', false)

    // Flag logic mirrors bamCoverage --filterRNAstrand:
    //   forward-stranded (F/FR): plus = flag16 NOT set (-F 2064); minus = flag16 set (-f 16 -F 2048)
    //   reverse-stranded:        plus = flag16 set    (-f 16 -F 2048); minus = flag16 NOT set (-F 2064)
    script:
    """
    if [ "${strand_arg}" == "F" ] || [ "${strand_arg}" == "FR" ]; then
        samtools view -@ ${task.cpus} -b -F 2064        -o ${sample}.ribo.plus.bam  ${bam}
        samtools view -@ ${task.cpus} -b -f 16 -F 2048  -o ${sample}.ribo.minus.bam ${bam}
    else
        samtools view -@ ${task.cpus} -b -f 16 -F 2048  -o ${sample}.ribo.plus.bam  ${bam}
        samtools view -@ ${task.cpus} -b -F 2064        -o ${sample}.ribo.minus.bam ${bam}
    fi
    samtools index -@ ${task.cpus} ${sample}.ribo.plus.bam
    samtools index -@ ${task.cpus} ${sample}.ribo.minus.bam
    """
}

///////////////////////////////////////////////////////////////////////////////
//  E N D   R I B O - S E Q   B I G W I G   G E N E R A T I O N
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  E N D   M E R G E D   S A M P L E   D E D U P L I C A T I O N   S T A T S
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//          E N D   G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////

// END OF GENOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

CLIP_LOG.map { sample, index, clip_log -> [ [sample, index], clip_log ] }
        .set { CLIP_LOG_INDEXED_FOR_GENOME }
FILTER_LOG.map { sample, index, filter_log -> [ [sample, index], filter_log ] }
          .set { FILTER_LOG_INDEXED_FOR_GENOME }

///////////////////////////////////////////////////////////////////////////////
// GENOME STATS

GENOME_ALIGNMENT_LOG_STATS
  .map { sample, index, log_file -> [ [sample, index], log_file] }
  .set { GENOME_ALIGNMENT_LOG_STATS_INDEXED }

GENOME_ALIGNMENT_SECONDARY_COUNT
  .map { sample, index, secondary -> [ [sample, index], secondary ] }
  .set { GENOME_ALIGNMENT_SECONDARY_COUNT_INDEXED }

// QPASS and per-lane dedup counts are 3-file tuples: total/primary/secondary.
// The shape change propagates into the join below and into the heredoc.
GENOME_QPASS_COUNTS
  .map { sample, index, total, primary, secondary ->
      [ [sample, index], total, primary, secondary ] }
  .set { GENOME_QPASS_COUNTS_INDEXED }

GENOME_INDIVIDUAL_DEDUP_COUNT
  .map { sample, index, total, primary, secondary ->
      [ [sample, index], total, primary, secondary ] }
  .set { GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED }

CLIP_LOG_INDEXED_FOR_GENOME.join(FILTER_LOG_INDEXED_FOR_GENOME)
            .join(GENOME_ALIGNMENT_LOG_STATS_INDEXED)
            .join(GENOME_ALIGNMENT_SECONDARY_COUNT_INDEXED)
            .join(GENOME_QPASS_COUNTS_INDEXED)
            .join(GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED)
            .map { key, clip_log, filter_log, genome_log, genome_sec,
                   qpass_total, qpass_primary, qpass_secondary,
                   dedup_total, dedup_primary, dedup_secondary ->
                [key[0], key[1], clip_log, filter_log, genome_log, genome_sec,
                 qpass_total, qpass_primary, qpass_secondary,
                 dedup_total, dedup_primary, dedup_secondary]
            }
            .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }

process individual_genome_alignment_stats {
    storeDir get_storedir('stats') + '/genome/individual'

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

// Stage dedup_* under fixed names so they never collide with qpass_* when
// dedup_method == 'none' aliases the qpass count channel into the dedup slot
// (the two channels would otherwise emit the same `*.qpass.*.count` files).
input:
    set val(sample), val(index),
        file(clip_log), file(filter_log), file(genome_log), file(genome_secondary_count),
        file(qpass_total), file(qpass_primary), file(qpass_secondary),
        file('dedup.total.count'), file('dedup.primary.count'), file('dedup.secondary.count') \
        from GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT

output:
    set val(sample), val(index), file("${sample}.${index}.genome_individual.csv") \
   into GENOME_INDIVIDUAL_ALIGNMENT_STATS

    """
    python3 - << 'PYEOF'
# cutadapt log
total_reads = 0; clipped_reads = 0
with open('${clip_log}') as fh:
    for line in fh:
        if line.startswith('Total reads'):
            total_reads = int(''.join(line.split()[-1].split(',')))
        elif line.startswith('Reads written'):
            clipped_reads = int(''.join(line.split()[-2].split(',')))

# filter log (bowtie2)
fl = [l for l in open('${filter_log}') if l.strip() and not l[0].isalpha()]
filtered_out = int(fl[3].split()[0]) + int(fl[4].split()[0])
filter_kept  = int(fl[2].split()[0])

# STAR genome alignment log — parsed by rfc parse-star-log
import csv, subprocess
proc = subprocess.run(
    ['rfc', 'parse-star-log', '${genome_log}'],
    check=True, capture_output=True, text=True,
)
star_row = next(csv.DictReader(proc.stdout.splitlines(), delimiter='\\t'))
genome_once  = int(star_row['uniquely_mapped'])
genome_multi = int(star_row['multi_loci_mapped'])
genome_unal  = int(star_row['unmapped_total'])
genome_primary = genome_once + genome_multi  # one primary record per mapped read

def read_int(p):
    with open(p) as fh:
        return int(fh.read().strip().split()[0])

genome_secondary = read_int('${genome_secondary_count}')
genome_total     = genome_primary + genome_secondary

qpass_total_v     = read_int('${qpass_total}')
qpass_primary_v   = read_int('${qpass_primary}')
qpass_secondary_v = read_int('${qpass_secondary}')

dedup_total_v     = read_int('dedup.total.count')
dedup_primary_v   = read_int('dedup.primary.count')
dedup_secondary_v = read_int('dedup.secondary.count')

# Raw count rows only — percentage rows are inserted at the combine stage by
# scripts/stats_percentage.py, because summing per-lane percentages would
# produce incorrect merged values.
rows = [
    ('total_reads',                  total_reads),
    ('clipped_reads',                clipped_reads),
    ('filtered_out',                 filtered_out),
    ('filter_kept',                  filter_kept),
    ('genome_aligned_once',          genome_once),
    ('genome_aligned_many',          genome_multi),
    ('genome_unaligned',             genome_unal),
    ('genome_primary_alignments',    genome_primary),
    ('genome_secondary_alignments',  genome_secondary),
    ('genome_total_alignments',      genome_total),
    ('qpass_primary_alignments',     qpass_primary_v),
    ('qpass_secondary_alignments',   qpass_secondary_v),
    ('qpass_total_alignments',       qpass_total_v),
    ('dedup_primary_alignments',     dedup_primary_v),
    ('dedup_secondary_alignments',   dedup_secondary_v),
    ('dedup_total_alignments',       dedup_total_v),
]
with open('${sample}.${index}.genome_individual.csv', 'w') as fh:
    fh.write(',${sample}.${index}\\n')
    for k, v in rows:
        fh.write(f'{k},{v}\\n')
PYEOF
    """
}

GENOME_INDIVIDUAL_ALIGNMENT_STATS
  .into { GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
      GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GROUPING}

GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
  .map { sample, index, stats_file -> stats_file }
  .toSortedList().set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED }

process combine_individual_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/individual'

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

input:
    file(stat_table) from GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

output:
    file('genome_individual_essential.csv') \
      into COMBINED_INDIVIDUAL_GENOME_ALIGNMENT_STATS

    script:
    if (stat_table.size() == 0) {
        '''
        echo "No individual statistics data available" > genome_individual_essential.csv
        '''
    } else {
        """
        rfc merge overall-stats \\
            -o raw_combined_individual_genome_aln_stats.csv ${stat_table}
        rfc genome-stats-percentage \\
            -i raw_combined_individual_genome_aln_stats.csv \\
            -o genome_individual_essential.csv
        """
    }
}

GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GROUPING
  .map { sample, index, file -> [ sample, file ] }
  .groupTuple()
  .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED }

// Merged dedup counts come from the dedup step (position or umicollapse). For
// dedup_method=='none' the channel is empty and the helper script falls back
// to summed per-lane qpass counts (which already populate the dedup rows).
if (dedup_method == 'position' || dedup_method == 'umicollapse') {
    GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
        .combine(GENOME_MERGED_DEDUP_COUNTS, by: 0)
        .map { sample, stat_files, total, primary, secondary ->
            [sample, stat_files, total, primary, secondary]
        }
        .set { GENOME_STATS_INPUT }
} else {
    GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
        .map { sample, stat_files -> [sample, stat_files, null, null, null] }
        .set { GENOME_STATS_INPUT }
}

process sum_individual_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/merged'

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

input:
    set val(sample), file(stat_files),
        val(merged_dedup_total),
        val(merged_dedup_primary),
        val(merged_dedup_secondary) from GENOME_STATS_INPUT

output:
    set val(sample), file("${sample}.genome_merged.csv") into GENOME_MERGED_ALIGNMENT_STATS

    script:
    has_merged_counts = merged_dedup_total != null && merged_dedup_total.toString() != 'null'

    if ((dedup_method == 'position' || dedup_method == 'umicollapse') && has_merged_counts) {
        """
        rfc sum-stats -n ${sample} -o ${sample}.genome_merged.tmp.csv ${stat_files}
        rfc update-dedup-counts \\
            --dedup-total-file ${merged_dedup_total} \\
            --dedup-primary-file ${merged_dedup_primary} \\
            --dedup-secondary-file ${merged_dedup_secondary} \\
            --input-csv ${sample}.genome_merged.tmp.csv \\
            --output-csv ${sample}.genome_merged.csv
        """
    } else {
        """
        rfc sum-stats -n ${sample} -o ${sample}.genome_merged.csv ${stat_files}
        """
    }
}

// Collect merged genome stats
GENOME_MERGED_ALIGNMENT_STATS
  .map { sample, stats_file -> stats_file }
  .toSortedList().set { GENOME_MERGED_ALIGNMENT_STATS_COLLECTED }

// COMBINE MERGED GENOME ALIGNMENT STATS
process combine_merged_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/merged'

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

input:
    file(stat_files) from GENOME_MERGED_ALIGNMENT_STATS_COLLECTED

output:
    file('genome_merged_essential.csv') into COMBINED_MERGED_GENOME_ALIGNMENT_STATS

    script:
    if (stat_files.size() == 0) {
        '''
        echo "No merged genome statistics data available" > genome_merged_essential.csv
        '''
    } else {
        """
        rfc merge overall-stats \\
            -o raw_combined_merged_genome_aln_stats.csv ${stat_files}
        rfc genome-stats-percentage \\
            -i raw_combined_merged_genome_aln_stats.csv \\
            -o genome_merged_essential.csv
        """
    }
}

// END OF GENOME STATS
///////////////////////////////////////////////////////////////////////////////


// Genome is now the only alignment method - use genome stats as final
COMBINED_INDIVIDUAL_GENOME_ALIGNMENT_STATS.set { FINAL_INDIVIDUAL_STATS }
COMBINED_MERGED_GENOME_ALIGNMENT_STATS.set { FINAL_MERGED_STATS }

// Use final stats directly for publishing
FINAL_INDIVIDUAL_STATS.set { FINAL_INDIVIDUAL_STATS_FOR_PUBLISH }
FINAL_MERGED_STATS.set { FINAL_MERGED_STATS_FOR_PUBLISH }

process publish_stats {
    publishDir get_publishdir('stats'), mode: 'copy'

    executor 'local'

  input:
    file(individual_stats) from FINAL_INDIVIDUAL_STATS_FOR_PUBLISH
    file(merged_stats) from FINAL_MERGED_STATS_FOR_PUBLISH

  output:
    file('individual_stats.csv') into INDIVIDUAL_STATS_PUBLISHED
    file('stats.csv') into MERGED_STATS_PUBLISHED

    """
  cp ${individual_stats} individual_stats.csv
  cp ${merged_stats} stats.csv
  """
}

////////////////////////////////////////////////////////////////////////////////
///////                       /* RNA-Seq */                            /////////
////////////////////////////////////////////////////////////////////////////////

do_rnaseq = params.get('do_rnaseq', false) && \
            params.get('rnaseq', false)

// This outer if clause contains the rest of the RNASEQ
if (do_rnaseq) {
    rnaseq_fastq_base =  params.get('rnaseq', [:]).get('fastq_base', '')
    if (! rnaseq_fastq_base.endsWith('/') && rnaseq_fastq_base != '') {
        rnaseq_fastq_base = "${rnaseq_fastq_base}/"
    }

    // RNA-seq accepts BOTH single-end and paired-end lanes. YAML schema per lane:
    //   - string                    -> single-end (one FASTQ)
    //   - [R1, R2] two-element list -> paired-end (R1/R2 mates with matching IDs)
    //
    // Channel tuple shape: [sample, index, reads] where `reads` is a list of
    // staged files: [r1] for SE, [r1, r2] for PE. Downstream processes branch
    // on `reads.size()` to assemble the right tool flags. This is the same
    // pattern the upstream nf-core modules use for read-input variability.
    def rnaseq_lane_tuples = []
    params.get('rnaseq', [:]).fastq.each { sample_name, lanes ->
        lanes.eachWithIndex { entry, i ->
            def lane_idx = i + 1
            def lane_files
            if ((entry instanceof List) && entry.size() == 2) {
                lane_files = [
                    file("${rnaseq_fastq_base}${entry[0]}"),
                    file("${rnaseq_fastq_base}${entry[1]}")]
            } else if (entry instanceof CharSequence) {
                lane_files = [file("${rnaseq_fastq_base}${entry}")]
            } else {
                throw new IllegalArgumentException(
                    "rnaseq.fastq.${sample_name} entry #${lane_idx} must be either:\n" +
                    "  * a string (single-end FASTQ path), or\n" +
                    "  * a 2-element list [R1, R2] (paired-end mates with matching read IDs).\n" +
                    "Got: ${entry}")
            }
            rnaseq_lane_tuples << [sample_name, lane_idx, lane_files]
        }
    }

    Channel.from(rnaseq_lane_tuples)
    .into {  RNASEQ_FASTQ_FASTQC
             RNASEQ_FASTQ_CLIP
             RNASEQ_FASTQ_EXISTENCE }

    if (params.do_check_file_existence) {
        RNASEQ_FASTQ_EXISTENCE
            .map { sample, index, reads ->
                reads.every { file_exists(it) }
            }
    }

    process rnaseq_raw_fastqc {
        publishDir get_publishdir('fastqc', true), mode: 'copy'

      input:
        set val(sample), val(index), file(reads) from RNASEQ_FASTQ_FASTQC

      output:
        file("*_fastqc.html") into RNASEQ_FASTQC_HTML
        file("*_fastqc.zip")  into RNASEQ_FASTQC_OUT

    when:
    params.do_fastqc && do_rnaseq

        // SE: `reads` is [r1]; PE: `reads` is [r1, r2]. fastqc takes any
        // number of input fastqs, so we just pass the whole list.
        script:
        """
        fastqc ${reads} --outdir=\$PWD -t ${task.cpus}
        """
    }

    process rnaseq_clip {
        storeDir get_storedir('clip', true)

  input:
        set val(sample), val(index), file(reads) from RNASEQ_FASTQ_CLIP

  output:
        // Glob captures 1 file for SE (R1 only) or 2 files for PE (R1+R2).
        // Downstream `file(reads)` will receive a list of the same shape.
        set val(sample), val(index),
            file("${sample}.${index}.clipped_R*.fastq.gz") into RNASEQ_CLIP_OUT
        set val(sample), val(index), file("${sample}.${index}.clipped.log") \
                                                      into RNASEQ_CLIP_LOG

        // SE: cutadapt with just -o; PE: cutadapt with -o + -p. The clipping
        // arguments (params.rnaseq.clip_arguments) apply to both modes.
        script:
        if (reads.size() == 2) {
            """
            cutadapt --cores=${task.cpus} ${params.rnaseq.clip_arguments} \\
                -o ${sample}.${index}.clipped_R1.fastq.gz \\
                -p ${sample}.${index}.clipped_R2.fastq.gz \\
                ${reads[0]} ${reads[1]} \\
                > ${sample}.${index}.clipped.log 2>&1
            """
        } else {
            """
            cutadapt --cores=${task.cpus} ${params.rnaseq.clip_arguments} \\
                -o ${sample}.${index}.clipped_R1.fastq.gz \\
                ${reads[0]} \\
                > ${sample}.${index}.clipped.log 2>&1
            """
        }
    }

    // Set RNA-seq clip log for genome stats
    RNASEQ_CLIP_LOG.set { RNASEQ_CLIP_LOG_FOR_GENOME }

    RNASEQ_FILTER_INDEX = Channel.from([[
             params.input.reference.filter
                .split('/')[-1]
                .replaceAll('\\*$', '')
                .replaceAll('\\.$', ''),
             file(params.input.reference.filter),
            ]])

    process rnaseq_filter {
        storeDir get_storedir('filter', true)

    input:
        set val(sample), val(index), file(reads) from RNASEQ_CLIP_OUT
        set val(bowtie2_index_base), file(bowtie2_index_files) \
                          from RNASEQ_FILTER_INDEX.first()

    output:
        set val(sample), val(index), file("${sample}.${index}.filter.bam") \
        into RNASEQ_FILTER_BAM
        set val(sample), val(index),
            file("${sample}.${index}.aligned.filter_R*.fastq.gz") into RNASEQ_FILTER_ALIGNED
        set val(sample), val(index),
            file("${sample}.${index}.unaligned.filter_R*.fastq.gz") into RNASEQ_FILTER_UNALIGNED
        set val(sample), val(index), file("${sample}.${index}.filter.log") \
        into RNASEQ_FILTER_LOG

        // bowtie2 SAM piped directly into samtools sort (no intermediate
        // `samtools view -b -`). The unused idxstats sidecar is gone.
        //
        // SE uses --al-gz / --un-gz with a single output path.
        // PE uses --al-conc-gz / --un-conc-gz with a `%` template that bowtie2
        // expands to 1/2 for the two mates.
        //
        // Decouple alignment threads from sort threads: bowtie2 doesn't scale
        // linearly past ~16 threads on the rRNA filter index, and samtools
        // sort plateaus around 8 threads. Capping both keeps the worst-case
        // heap footprint within task.memory and avoids OOM-killing bowtie2
        // (which previously segfaulted with exit 139 under memory pressure).
        script:
        def aln_threads  = Math.min(task.cpus as int, 16)
        def sort_threads = Math.min(task.cpus as int, 8)
        def sort_mem     = samtools_sort_mem_per_thread_mb(task)
        if (reads.size() == 2) {
            """
            set -o pipefail
            bowtie2 ${params.rnaseq.filter_arguments} \\
                    -x ${bowtie2_index_base} -1 ${reads[0]} -2 ${reads[1]} \\
                    --threads ${aln_threads} \\
                    --al-conc-gz ${sample}.${index}.aligned.filter_R%.fastq.gz \\
                    --un-conc-gz ${sample}.${index}.unaligned.filter_R%.fastq.gz \\
                    2> ${sample}.${index}.filter.log \\
                | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${sample}.${index}.filter.bam -
            samtools index -@ ${sort_threads} ${sample}.${index}.filter.bam
            """
        } else {
            """
            set -o pipefail
            bowtie2 ${params.rnaseq.filter_arguments} \\
                    -x ${bowtie2_index_base} -U ${reads[0]} \\
                    --threads ${aln_threads} \\
                    --al-gz ${sample}.${index}.aligned.filter_R1.fastq.gz \\
                    --un-gz ${sample}.${index}.unaligned.filter_R1.fastq.gz \\
                    2> ${sample}.${index}.filter.log \\
                | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${sample}.${index}.filter.bam -
            samtools index -@ ${sort_threads} ${sample}.${index}.filter.bam
            """
        }
    }

    RNASEQ_FILTER_LOG.set { RNASEQ_FILTER_LOG_FOR_GENOME }
    RNASEQ_FILTER_UNALIGNED.set { RNASEQ_FILTER_UNALIGNED_GENOME }

///////////////////////////////////////////////////////////////////////////////
//          R N A - S E Q   G E N O M E   A L I G N M E N T
///////////////////////////////////////////////////////////////////////////////

    // RNA-seq genome alignment — feed rRNA-filtered reads directly to STAR
    RNASEQ_GENOME_INPUT_CHANNEL_RAW = RNASEQ_FILTER_UNALIGNED_GENOME

    RNASEQ_GENOME_INDEX = Channel.from([[ "star_index", file(params.input.reference.genome) ]])
    RNASEQ_GENOME_INDEX.set { RNASEQ_GENOME_INDEX_FOR_ALIGN }

    RNASEQ_GENOME_INPUT_CHANNEL_RAW.set { RNASEQ_GENOME_INPUT_CHANNEL }

    process rnaseq_genome_alignment {
        storeDir get_storedir('genome_alignment', true) + '/' + params.output.individual_lane_directory

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

    input:
        set val(sample), val(index), file(reads) from RNASEQ_GENOME_INPUT_CHANNEL
        set val(genome_base), file(genome_files) from RNASEQ_GENOME_INDEX_FOR_ALIGN.first()

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.bam") \
        into RNASEQ_GENOME_ALIGNMENT_BAM
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.log") \
        into RNASEQ_GENOME_ALIGNMENT_LOG
        set val(sample), val(index),
            file("${sample}.${index}.rnaseq_genome_alignment.secondary.count") \
        into RNASEQ_GENOME_ALIGNMENT_SECONDARY_COUNT

        // Mirrors ENCODE rna-seq-pipeline alignment (align.py).
        // --outSAMunmapped Within keeps unmapped reads in the BAM.
        // STAR's --readFilesIn takes 1 path for SE or 2 space-separated paths
        // for PE — we build the right argument from the `reads` list.
        script:
        def read_files_in = reads.join(' ')
        """
    set -o pipefail
    mkdir -p star_out

    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${read_files_in} \\
        --readFilesCommand zcat \\
        --readFilesType Fastx \\
        --genomeDir ${genome_files} \\
        --genomeLoad NoSharedMemory \\
        ${params.star.rnaseq_arguments} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes All \\
        --outSAMstrandField intronMotif \\
        --outSAMunmapped Within \\
        --outSAMattrRGline ID:${sample}.${index} SM:${sample} PL:ILLUMINA \\
        --outFileNamePrefix star_out/

    mv star_out/Aligned.sortedByCoord.out.bam ${sample}.${index}.rnaseq_genome_alignment.bam
    samtools index -@ ${task.cpus} ${sample}.${index}.rnaseq_genome_alignment.bam
    cp star_out/Log.final.out ${sample}.${index}.rnaseq_genome_alignment.log

    # Secondary count for the stats alignment breakdown.
    samtools view -@ ${task.cpus} -c -f 256 ${sample}.${index}.rnaseq_genome_alignment.bam \\
        > ${sample}.${index}.rnaseq_genome_alignment.secondary.count
        """
    }

    RNASEQ_GENOME_ALIGNMENT_BAM.set { RNASEQ_GENOME_ALIGNMENT_BAM_FOR_QUALITY }
    RNASEQ_GENOME_ALIGNMENT_LOG.set { RNASEQ_GENOME_ALIGNMENT_LOG_FOR_STATS }

    // RNA-seq genome quality filtering
    process rnaseq_genome_quality_filter {
        storeDir get_storedir('quality_filter', true)

    input:
        set val(sample), val(index), file(bam) from RNASEQ_GENOME_ALIGNMENT_BAM_FOR_QUALITY

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.qpass.bam") \
            into RNASEQ_GENOME_ALIGNMENT_QPASS_BAM
        set val(sample), val(index),
            file("${sample}.${index}.qpass.total.count"),
            file("${sample}.${index}.qpass.primary.count"),
            file("${sample}.${index}.qpass.secondary.count") \
            into RNASEQ_GENOME_QPASS_COUNTS_OUT

        // -F rnaseq.filter_flags: drop unmapped (4) + supplementary (2048);
        // by default RNA-seq keeps secondary records (256) so the alignment
        // breakdown (total/primary/secondary) is meaningful.
        script:
        def rnaseq_cfg   = params.get('rnaseq', [:])
        def rna_mapq     = rnaseq_cfg.get('mapping_quality_cutoff', 4)
        def rna_flags    = rnaseq_cfg.get('filter_flags', 2308)
        """
        samtools view -@ ${task.cpus} -bq ${rna_mapq} -F ${rna_flags} ${bam} \\
            > ${sample}.${index}.rnaseq_genome_alignment.qpass.bam
        samtools index -@ ${task.cpus} ${sample}.${index}.rnaseq_genome_alignment.qpass.bam
        samtools view -@ ${task.cpus} -c        ${sample}.${index}.rnaseq_genome_alignment.qpass.bam \\
            > ${sample}.${index}.qpass.total.count
        samtools view -@ ${task.cpus} -c -F 2304 ${sample}.${index}.rnaseq_genome_alignment.qpass.bam \\
            > ${sample}.${index}.qpass.primary.count
        samtools view -@ ${task.cpus} -c -f 256  ${sample}.${index}.rnaseq_genome_alignment.qpass.bam \\
            > ${sample}.${index}.qpass.secondary.count
        """
    }

    RNASEQ_GENOME_ALIGNMENT_QPASS_BAM.into {
        RNASEQ_GENOME_QPASS_BAM_FOR_BED
        RNASEQ_GENOME_QPASS_BAM_FOR_NODEDUP_MERGE
        RNASEQ_GENOME_QPASS_BAM_FOR_PUBLISH_NONE  // per-lane qpass BAM publish on 'none' path
    }

    // QPASS count tuples (total, primary, secondary) feed two sinks:
    //   1) the per-lane stats join (always)
    //   2) the dedup-count alias when rnaseq_dedup_method == 'none'
    // NF DSL1 forbids consuming the same channel twice, so we split here.
    RNASEQ_GENOME_QPASS_COUNTS_OUT.into {
        RNASEQ_GENOME_QPASS_COUNTS
        RNASEQ_GENOME_QPASS_COUNTS_FOR_NONE_DEDUP
    }

    // Convert individual RNA-seq qpass BAM → per-lane BED. Published as the
    // final per-lane artifact when rnaseq.dedup_method == 'none'.
    process rnaseq_individual_genome_bam_to_bed {
        storeDir get_storedir('bam_to_bed', true) + '/' + params.output.individual_lane_directory
        publishDir get_publishdir('alignments') + '/rnaseq/' + params.output.individual_lane_directory,
                   mode: 'copy', enabled: rnaseq_dedup_method == 'none'

    input:
        set val(sample), val(index), file(bam) from RNASEQ_GENOME_QPASS_BAM_FOR_BED

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.qpass.bed") \
            into RNASEQ_INDIVIDUAL_GENOME_BED

        """
        if [ \$(samtools view -@ ${task.cpus} -c ${bam}) -eq 0 ]; then
            touch ${sample}.${index}.rnaseq_genome.qpass.bed
        else
            bamToBed -i ${bam} > ${sample}.${index}.rnaseq_genome.qpass.bed
        fi
        """
    }

    RNASEQ_INDIVIDUAL_GENOME_BED.into {
        RNASEQ_GENOME_BED_FOR_INDEX_COL          // position dedup: tag with sample.index
        RNASEQ_GENOME_BED_FOR_INDEX_SEP_PRE      // position dedup: lane lookup for separation
        RNASEQ_GENOME_BED_FOR_NONE_MERGE         // dedup=none: concat to merged qpass BED
    }

    // Add sample index column (column 7) so the merged BED can be awk-split per
    // lane after rfc dedup. Position-dedup only.
    process rnaseq_add_sample_index_col_to_genome_bed {
        storeDir get_storedir('bam_to_bed', true) + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bed) from RNASEQ_GENOME_BED_FOR_INDEX_COL

    output:
        set val(sample), file("${sample}.${index}.rnaseq_genome.with_sample_index.bed") \
            into RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED

    when:
        rnaseq_dedup_method == 'position'

        """
        awk -v newcol=${sample}.${index} '{print(\$0"\\t"newcol)}' ${bed} \\
            > ${sample}.${index}.rnaseq_genome.with_sample_index.bed
        """
    }

    RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED.groupTuple()
    .set { RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED }

    process rnaseq_merge_genome_bed {
        storeDir get_storedir('bam_to_bed', true) + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(bed_files) from RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED

    output:
        set val(sample), file("${sample}.merged.pre_dedup.bed") \
            into RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP

    when:
        rnaseq_dedup_method == 'position'

        """
        cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.merged.pre_dedup.bed
        """
    }

    // Index per-step inputs (clip log, filter log, genome log, secondary count,
    // qpass count tuple) by (sample, lane) for the per-lane stats join.
    RNASEQ_CLIP_LOG_FOR_GENOME.map { sample, index, clip_log -> [ [sample, index], clip_log ] }
        .set { RNASEQ_CLIP_LOG_INDEXED_FOR_GENOME }
    RNASEQ_FILTER_LOG_FOR_GENOME.map { sample, index, filter_log -> [ [sample, index], filter_log ] }
        .set { RNASEQ_FILTER_LOG_INDEXED_FOR_GENOME }
    RNASEQ_GENOME_ALIGNMENT_LOG_FOR_STATS.map { sample, index, genome_log -> [ [sample, index], genome_log ] }
        .set { RNASEQ_GENOME_ALIGNMENT_LOG_INDEXED }
    RNASEQ_GENOME_ALIGNMENT_SECONDARY_COUNT.map { sample, index, secondary -> [ [sample, index], secondary ] }
        .set { RNASEQ_GENOME_ALIGNMENT_SECONDARY_COUNT_INDEXED }
    RNASEQ_GENOME_QPASS_COUNTS.map { sample, index, total, primary, secondary ->
            [ [sample, index], total, primary, secondary ] }
        .set { RNASEQ_GENOME_QPASS_COUNTS_INDEXED }

///////////////////////////////////////////////////////////////////////////////
//          R N A - S E Q   G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////
//
// rnaseq_merge_genome_alignment (unfiltered merge) and
// rnaseq_genome_create_nodedup_counts (duplicate wc -l) have been removed —
// neither had downstream consumers.

    // Merge quality-filtered per-lane BAMs to a single merged qpass BAM. This
    // is the bigwig source for dedup_method == 'none' AND is published in that
    // case. For position dedup it stays in intermediates and feeds the
    // post-dedup BAM extraction.
    RNASEQ_GENOME_QPASS_BAM_FOR_NODEDUP_MERGE
        .map { sample, index, bam -> [sample, bam] }
        .groupTuple()
        .set { RNASEQ_GENOME_QPASS_GROUPED_BAM }

    process rnaseq_merge_quality_filtered_bam {
        storeDir get_storedir('genome_alignment', true) + '/' + params.output.merged_lane_directory
        publishDir get_publishdir('alignments') + '/rnaseq/' + params.output.merged_lane_directory,
                   mode: 'copy', enabled: rnaseq_dedup_method == 'none'

        input:
        set val(sample), file(bam_files) from RNASEQ_GENOME_QPASS_GROUPED_BAM

        output:
        set val(sample),
            file("${sample}.rnaseq_genome.qpass.merged.bam"),
            file("${sample}.rnaseq_genome.qpass.merged.bam.bai") \
            into RNASEQ_GENOME_QPASS_MERGED_BAM

        """
        samtools merge -@ ${task.cpus} ${sample}.rnaseq_genome.qpass.merged.bam ${bam_files}
        samtools index -@ ${task.cpus} ${sample}.rnaseq_genome.qpass.merged.bam
        """
    }

    RNASEQ_GENOME_QPASS_MERGED_BAM.into {
        RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_NODEDUP
        RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_DEDUP
    }

    // RNA-seq genome deduplication (position-based only). When
    // rnaseq.dedup_method == 'none' the position chain is skipped via `when:`.
    process rnaseq_genome_deduplicate {
        storeDir get_publishdir('alignments') + '/rnaseq/' + params.output.merged_lane_directory

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

        input:
        set val(sample), file(bed) from RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP

        output:
        set val(sample), file("${sample}.rnaseq_genome.post_dedup.bed") \
            into RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP

        when:
        rnaseq_dedup_method == 'position'

        """
        sort -k1,1 -k2,2n -k3,3n ${bed} > ${sample}.rnaseq_genome.sorted.bed
        rfc dedup -i ${sample}.rnaseq_genome.sorted.bed -o ${sample}.rnaseq_genome.post_dedup.bed
        """
    }

    if (rnaseq_dedup_method == 'position') {

        // Split the merged dedup BED for per-lane separation and for the
        // sample-level BAM extraction (which becomes the bigwig source).
        RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.into {
            RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_SEP
            RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_BIGWIG
        }

        RNASEQ_GENOME_BED_FOR_INDEX_SEP_PRE
            .map { sample, index, file -> [sample, index] }
            .combine(RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_SEP, by: 0)
            .set { RNASEQ_GENOME_BED_FOR_INDEX_SEP_POST_DEDUP }

        // Awk-split per lane and emit total/primary/secondary count tuple
        // derived from the per-lane BED (total = lines, primary = unique read
        // names, secondary = total - primary).
        process rnaseq_separate_genome_bed_post_dedup {
            storeDir get_publishdir('alignments') + '/rnaseq/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(bed) from RNASEQ_GENOME_BED_FOR_INDEX_SEP_POST_DEDUP

        output:
            set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.post_dedup.bed") \
                into RNASEQ_GENOME_BED_DEDUPLICATED
            set val(sample), val(index),
                file("${sample}.${index}.dedup.total.count"),
                file("${sample}.${index}.dedup.primary.count"),
                file("${sample}.${index}.dedup.secondary.count") \
                into RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_POSITION

        when:
            rnaseq_dedup_method == 'position'

            """
            awk -v s=${sample}.${index} \\
                '\$7 == s { print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6 }' \\
                ${bed} > ${sample}.${index}.rnaseq_genome.post_dedup.bed
            total=\$(wc -l < ${sample}.${index}.rnaseq_genome.post_dedup.bed)
            primary=\$(awk '{print \$4}' ${sample}.${index}.rnaseq_genome.post_dedup.bed | sort -u | wc -l)
            secondary=\$((total - primary))
            echo \${total}     > ${sample}.${index}.dedup.total.count
            echo \${primary}   > ${sample}.${index}.dedup.primary.count
            echo \${secondary} > ${sample}.${index}.dedup.secondary.count
            """
        }

        RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_POSITION
            .set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT }

        // Sample-level: merged BED → merged post-dedup BAM via read extraction
        // from the qpass merged BAM, then emit alignment-count breakdown from
        // the BAM. This bigwig source replaces the qpass merged BAM in the
        // 'position' branch (see rnaseq_create_bigwig wiring below).
        RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_BIGWIG
            .join(RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_DEDUP)
            .map { sample, bed, qpass_bam, qpass_bai -> [sample, bed, qpass_bam] }
            .set { RNASEQ_GENOME_DEDUP_BED_WITH_QPASS_BAM }

        process rnaseq_extract_dedup_reads_from_bam {
            storeDir get_publishdir('alignments') + '/rnaseq/' + params.output.merged_lane_directory

            input:
            set val(sample), file(post_dedup_bed), file(qpass_bam) from RNASEQ_GENOME_DEDUP_BED_WITH_QPASS_BAM

            output:
            set val(sample),
                file("${sample}.rnaseq_genome.post_dedup.bam"),
                file("${sample}.rnaseq_genome.post_dedup.bam.bai") \
                into RNASEQ_GENOME_DEDUP_BAM_FOR_BIGWIG
            set val(sample),
                file("${sample}.merged_dedup.total.count"),
                file("${sample}.merged_dedup.primary.count"),
                file("${sample}.merged_dedup.secondary.count") \
                into RNASEQ_GENOME_MERGED_DEDUP_COUNTS

            when:
            rnaseq_dedup_method == 'position'

            """
            set +u
            source \$(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome
            set -u
            rfc extract-dedup-reads \\
                --bam ${qpass_bam} \\
                --bed ${post_dedup_bed} \\
                --output ${sample}.rnaseq_genome.post_dedup.bam
            samtools index -@ ${task.cpus} ${sample}.rnaseq_genome.post_dedup.bam
            samtools view -@ ${task.cpus} -c        ${sample}.rnaseq_genome.post_dedup.bam \\
                > ${sample}.merged_dedup.total.count
            samtools view -@ ${task.cpus} -c -F 2304 ${sample}.rnaseq_genome.post_dedup.bam \\
                > ${sample}.merged_dedup.primary.count
            samtools view -@ ${task.cpus} -c -f 256  ${sample}.rnaseq_genome.post_dedup.bam \\
                > ${sample}.merged_dedup.secondary.count
            """
        }

        RNASEQ_GENOME_DEDUP_BAM_FOR_BIGWIG.set { RNASEQ_GENOME_FOR_BIGWIG }

    } else {

        // rnaseq.dedup_method == 'none': dedup count tuple aliases qpass,
        // bigwig comes from the merged qpass BAM, and we publish a merged
        // qpass BED (concat of per-lane qpass BEDs).
        RNASEQ_GENOME_QPASS_COUNTS_FOR_NONE_DEDUP
            .set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT }
        RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_NODEDUP.set { RNASEQ_GENOME_FOR_BIGWIG }
        Channel.empty().set { RNASEQ_GENOME_MERGED_DEDUP_COUNTS }

        RNASEQ_GENOME_BED_FOR_NONE_MERGE
            .map { sample, index, bed -> [sample, bed] }
            .groupTuple()
            .set { RNASEQ_GENOME_BED_FOR_NONE_MERGE_GROUPED }

        process rnaseq_concat_genome_qpass_bed_for_publish_none {
            storeDir get_storedir('bam_to_bed', true) + '/' + params.output.merged_lane_directory
            publishDir get_publishdir('alignments') + '/rnaseq/' + params.output.merged_lane_directory,
                       mode: 'copy'

            input:
            set val(sample), file(bed_files) from RNASEQ_GENOME_BED_FOR_NONE_MERGE_GROUPED

            output:
            set val(sample), file("${sample}.rnaseq_genome.qpass.bed") \
                into RNASEQ_MERGED_QPASS_BED_NONE

            when:
            rnaseq_dedup_method == 'none'

            """
            cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.rnaseq_genome.qpass.bed
            """
        }
    }

    // Build the per-lane stats input.
    RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT
        .map { sample, index, total, primary, secondary ->
            [ [sample, index], total, primary, secondary ] }
        .set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED }

    RNASEQ_CLIP_LOG_INDEXED_FOR_GENOME.join(RNASEQ_FILTER_LOG_INDEXED_FOR_GENOME)
            .join(RNASEQ_GENOME_ALIGNMENT_LOG_INDEXED)
            .join(RNASEQ_GENOME_ALIGNMENT_SECONDARY_COUNT_INDEXED)
            .join(RNASEQ_GENOME_QPASS_COUNTS_INDEXED)
            .join(RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED)
            .map { key, clip_log, filter_log, genome_log, genome_sec,
                   qpass_total, qpass_primary, qpass_secondary,
                   dedup_total, dedup_primary, dedup_secondary ->
                [key[0], key[1], clip_log, filter_log, genome_log, genome_sec,
                 qpass_total, qpass_primary, qpass_secondary,
                 dedup_total, dedup_primary, dedup_secondary]
            }
            .set { RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }

    process rnaseq_individual_genome_alignment_stats {
        storeDir get_storedir('stats', true) + '/genome/individual'

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

        // Stage dedup_* under fixed names so they never collide with qpass_* when
        // rnaseq_dedup_method == 'none' aliases the qpass count channel into the
        // dedup slot (the two channels would otherwise emit the same files).
        input:
        set val(sample), val(index),
            file(clip_log), file(filter_log), file(genome_log), file(genome_secondary_count),
            file(qpass_total), file(qpass_primary), file(qpass_secondary),
            file('dedup.total.count'), file('dedup.primary.count'), file('dedup.secondary.count') \
            from RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT

        output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_individual.csv") \
            into RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS

        """
        python3 - << 'PYEOF'
# cutadapt log
total_reads = 0; clipped_reads = 0
with open('${clip_log}') as fh:
    for line in fh:
        if line.startswith('Total reads'):
            total_reads = int(''.join(line.split()[-1].split(',')))
        elif line.startswith('Reads written'):
            clipped_reads = int(''.join(line.split()[-2].split(',')))

# filter log (bowtie2) — paired-end, lines indexed differently than ribo SE.
fl = [l for l in open('${filter_log}') if l.strip() and not l[0].isalpha()]
filtered_out = int(fl[3].split()[0]) + int(fl[4].split()[0])
filter_kept  = int(fl[2].split()[0])

# STAR genome alignment log
import csv, subprocess
proc = subprocess.run(
    ['rfc', 'parse-star-log', '${genome_log}'],
    check=True, capture_output=True, text=True,
)
star_row = next(csv.DictReader(proc.stdout.splitlines(), delimiter='\\t'))
genome_once  = int(star_row['uniquely_mapped'])
genome_multi = int(star_row['multi_loci_mapped'])
genome_unal  = int(star_row['unmapped_total'])
genome_primary = genome_once + genome_multi

def read_int(p):
    with open(p) as fh:
        return int(fh.read().strip().split()[0])

genome_secondary = read_int('${genome_secondary_count}')
genome_total     = genome_primary + genome_secondary

qpass_total_v     = read_int('${qpass_total}')
qpass_primary_v   = read_int('${qpass_primary}')
qpass_secondary_v = read_int('${qpass_secondary}')

dedup_total_v     = read_int('dedup.total.count')
dedup_primary_v   = read_int('dedup.primary.count')
dedup_secondary_v = read_int('dedup.secondary.count')

rows = [
    ('total_reads',                  total_reads),
    ('clipped_reads',                clipped_reads),
    ('filtered_out',                 filtered_out),
    ('filter_kept',                  filter_kept),
    ('genome_aligned_once',          genome_once),
    ('genome_aligned_many',          genome_multi),
    ('genome_unaligned',             genome_unal),
    ('genome_primary_alignments',    genome_primary),
    ('genome_secondary_alignments',  genome_secondary),
    ('genome_total_alignments',      genome_total),
    ('qpass_primary_alignments',     qpass_primary_v),
    ('qpass_secondary_alignments',   qpass_secondary_v),
    ('qpass_total_alignments',       qpass_total_v),
    ('dedup_primary_alignments',     dedup_primary_v),
    ('dedup_secondary_alignments',   dedup_secondary_v),
    ('dedup_total_alignments',       dedup_total_v),
]
with open('${sample}.${index}.rnaseq_genome_individual.csv', 'w') as fh:
    fh.write(',${sample}.${index}\\n')
    for k, v in rows:
        fh.write(f'{k},{v}\\n')
PYEOF
        """
    }

///////////////////////////////////////////////////////////////////////////////
//  R N A - S E Q   B I G W I G   G E N E R A T I O N
///////////////////////////////////////////////////////////////////////////////
//
// Source switches by rnaseq.dedup_method:
//   - 'none'     → merged qpass BAM
//   - 'position' → merged post-dedup BAM (extracted via rfc-dedup BED)
// Both produce the same bigwig artifact name; the BAM source is the only diff.

    process rnaseq_create_bigwig {
        storeDir get_publishdir('bigwigs') + '/rnaseq'

        input:
        set val(sample), file(bam), file(bai) from RNASEQ_GENOME_FOR_BIGWIG

        output:
        set val(sample), file("${sample}.rnaseq.bigWig") into RNASEQ_BIGWIGS

        // bw_threads capped at 8 — see note on the ribo bigwig process. Each
        // `bamCoverage -p N` forks N Python workers; 32 per sample × multiple
        // samples easily exhausts the fork table and triggers fork() ENOMEM.
        script:
        def bw_threads = Math.min(task.cpus as int, 8)
        """
        bamCoverage -b ${bam} -o ${sample}.rnaseq.bigWig --outFileFormat bigwig \\
            --binSize 1 -p ${bw_threads}
        """
    }

///////////////////////////////////////////////////////////////////////////////
//  E N D   R N A - S E Q   B I G W I G   G E N E R A T I O N
///////////////////////////////////////////////////////////////////////////////

    // RNA-seq genome merged statistics compilation (following ribo-seq genome pattern)
    RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS.into {
        RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
        RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GROUPING
    }

    // Collect individual RNA-seq genome stats
    RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
    .map { sample, index, stats_file -> stats_file }
    .toSortedList().set { RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED }

    process combine_individual_rnaseq_genome_alignment_stats {
        executor 'local'
        storeDir get_storedir('stats', true) + '/genome/individual'
        publishDir get_rnaseq_publishdir('stats'), mode: 'copy'

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

    input:
        file(stat_table) from RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

    output:
        file('rnaseq_individual_stats.csv') \
          into COMBINED_INDIVIDUAL_RNASEQ_GENOME_ALIGNMENT_STATS

        script:
        if (stat_table.size() == 0) {
            '''
            echo "No RNA-seq individual statistics data available" > rnaseq_individual_stats.csv
            '''
        } else {
            """
            rfc merge overall-stats \\
                -o raw_combined_individual_rnaseq_genome_aln_stats.csv ${stat_table}
            rfc genome-stats-percentage \\
                -i raw_combined_individual_rnaseq_genome_aln_stats.csv \\
                -o rnaseq_individual_stats.csv
            """
        }
    }

    // Sum individual RNA-seq genome stats by sample
    RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GROUPING
    .map { sample, index, file -> [ sample, file ] }
    .groupTuple()
    .set { RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED }

    if (rnaseq_dedup_method == 'position') {
        RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
            .combine(RNASEQ_GENOME_MERGED_DEDUP_COUNTS, by: 0)
            .map { sample, stat_files, total, primary, secondary ->
                [sample, stat_files, total, primary, secondary] }
            .set { RNASEQ_GENOME_STATS_INPUT }
    } else {
        RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
            .map { sample, stat_files -> [sample, stat_files, null, null, null] }
            .set { RNASEQ_GENOME_STATS_INPUT }
    }

    process sum_individual_rnaseq_genome_alignment_stats {
        executor 'local'
        storeDir get_storedir('stats', true) + '/genome/merged'

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

        input:
        set val(sample), file(stat_files),
            val(merged_dedup_total),
            val(merged_dedup_primary),
            val(merged_dedup_secondary) from RNASEQ_GENOME_STATS_INPUT

        output:
        set val(sample), file("${sample}.rnaseq_genome_merged.csv") into RNASEQ_GENOME_MERGED_ALIGNMENT_STATS

        script:
        has_merged_counts = merged_dedup_total != null && merged_dedup_total.toString() != 'null'

        if (rnaseq_dedup_method == 'position' && has_merged_counts) {
            """
            rfc sum-stats -n ${sample} -o ${sample}.rnaseq_genome_merged.tmp.csv ${stat_files}
            rfc update-dedup-counts \\
                --dedup-total-file ${merged_dedup_total} \\
                --dedup-primary-file ${merged_dedup_primary} \\
                --dedup-secondary-file ${merged_dedup_secondary} \\
                --input-csv ${sample}.rnaseq_genome_merged.tmp.csv \\
                --output-csv ${sample}.rnaseq_genome_merged.csv
            """
        } else {
            """
            rfc sum-stats -n ${sample} -o ${sample}.rnaseq_genome_merged.csv ${stat_files}
            """
        }
    }

    // Combine merged RNA-seq genome stats
    RNASEQ_GENOME_MERGED_ALIGNMENT_STATS
    .map { sample, stats_file -> stats_file }
    .toSortedList().set { RNASEQ_GENOME_MERGED_ALIGNMENT_STATS_COLLECTED }

    process combine_merged_rnaseq_genome_alignment_stats {
        executor 'local'
        storeDir get_storedir('stats', true) + '/genome/merged'
        publishDir get_rnaseq_publishdir('stats'), mode: 'copy'

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

    input:
        file(stat_table) from RNASEQ_GENOME_MERGED_ALIGNMENT_STATS_COLLECTED

    output:
        file('rnaseq_stats.csv') \
          into COMBINED_MERGED_RNASEQ_GENOME_ALIGNMENT_STATS

        script:
        if (stat_table.size() == 0) {
            '''
            echo "No RNA-seq merged statistics data available" > rnaseq_stats.csv
            '''
        } else {
            """
            rfc merge overall-stats \\
                -o raw_combined_merged_rnaseq_genome_aln_stats.csv ${stat_table}
            rfc genome-stats-percentage \\
                -i raw_combined_merged_rnaseq_genome_aln_stats.csv \\
                -o rnaseq_stats.csv
            """
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    //          E N D   R N A - S E Q   G E N O M E   D E D U P L I C A T I O N
    ///////////////////////////////////////////////////////////////////////////////
    }
// RNA-Seq
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// (.ribo file creation has been removed — this is a pure genome-alignment pipeline.
//  Genome-aligned reads cannot be represented in the transcriptome-coordinate .ribo
//  format. Final artifacts are bigWigs and dedup BAM/BED files under output/.)
////////////////////////////////////////////////////////////////////////////////
