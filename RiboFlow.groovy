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

rnaseq_dedup_method = params.get('rnaseq', [:]).get('dedup_method', 'position').toString()

fastq_base = params.get('input', [:]).get('fastq_base', '')
if (! fastq_base.endsWith('/') && fastq_base != '') {
    fastq_base = "${fastq_base}/"
}

Channel.from(params.get('input', [:]).fastq.collect { k, v ->
    v.collect { z -> [k, v.indexOf(z) + 1,
                                                   file("${fastq_base}${z}")] }  })
    .flatten().collate(3).into {  INPUT_SAMPLES_VERBOSE
                             INPUT_SAMPLES_MD5
                                 INPUT_SAMPLES_EXISTENCE
                                 INPUT_SAMPLES_FASTQC
                             INPUT_SAMPLES_CLIP
                                 INPUT_SAMPLES_LOG
                             INPUT_SAMPLES_READ_LENGTH
                             INPUT_FOR_METADATA}

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

CLIP_OUT.into { CLIP_OUT_FASTQC; CLIP_OUT_DEDUP; CLIP_OUT_FILTER ; CLIP_OUT_READ_LENGTH }

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
    set val(sample), val(index), file("${sample}.${index}.filter.bam.bai") \
        into FILTER_BAI
    set val(sample), val(index), file("${sample}.${index}.aligned.filter.fastq.gz") \
        into FILTER_ALIGNED
    set val(sample), val(index), file("${sample}.${index}.unaligned.filter.fastq.gz") \
        into FILTER_UNALIGNED
    set val(sample), val(index), file("${sample}.${index}.filter.log") \
        into FILTER_LOG
    set val(sample), val(index), file("${sample}.${index}.filter.stats") \
        into FILTER_STATS

    """
    set -o pipefail
    bowtie2 ${params.alignment_arguments.filter} \
            -x ${bowtie2_index_base} -q ${fastq} \
            --threads ${task.cpus} \
            --al-gz ${sample}.${index}.aligned.filter.fastq.gz \
            --un-gz ${sample}.${index}.unaligned.filter.fastq.gz \
                     2> ${sample}.${index}.filter.log \
            | samtools view -@ ${task.cpus} -b - \
            | samtools sort -@ ${task.cpus} -o ${sample}.${index}.filter.bam \
            && samtools index -@ ${task.cpus} ${sample}.${index}.filter.bam \
            && samtools idxstats -@ ${task.cpus} ${sample}.${index}.filter.bam  > \
               ${sample}.${index}.filter.stats
    """
}

FILTER_ALIGNED.into { FILTER_ALIGNED_FASTQ_READ_LENGTH
    FILTER_ALIGNED_FASTQ_FASTQC}

FILTER_UNALIGNED.into { FILTER_UNALIGNED_FASTQ_READ_LENGTH
    FILTER_UNALIGNED_FASTQ_FASTQC
    FILTER_UNALIGNED_GENOME}

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
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.bam.bai") \
    into GENOME_ALIGNMENT_BAI
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.bam") optional true \
    into GENOME_TRANSCRIPTOME_BAM
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.aligned.fastq.gz") \
    into GENOME_ALIGNMENT_ALIGNED
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.unaligned.fastq.gz") \
    into GENOME_ALIGNMENT_UNALIGNED
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.log") \
    into GENOME_ALIGNMENT_LOG
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.stats") \
    into GENOME_ALIGNMENT_STATS
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.csv") \
    into GENOME_ALIGNMENT_CSV
    set val(sample), val(strand_arg) into GENOME_SAMPLE_STRAND

    script:
    def emit_tx_bam = params.get('star', [:]).get('output_transcriptome_bam', false)
    def quant_mode_arg = emit_tx_bam ? '--quantMode TranscriptomeSAM \\\n        ' : ''
    def tx_bam_cmd = emit_tx_bam \
        ? "mv star_out/Aligned.toTranscriptome.out.bam ${sample}.${index}.transcriptome_alignment.bam" \
        : ''
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
    samtools sort -@ ${task.cpus} \\
        -o ${sample}.${index}.genome_alignment.bam \\
        star_out/Aligned.sortedByCoord.out.bam
    rm -f star_out/Aligned.sortedByCoord.out.bam
    samtools index -@ ${task.cpus} ${sample}.${index}.genome_alignment.bam

    ${tx_bam_cmd}

    samtools fastq -@ ${task.cpus} -F 4 ${sample}.${index}.genome_alignment.bam \\
        | gzip > ${sample}.${index}.genome_alignment.aligned.fastq.gz

    if [ -f star_out/Unmapped.out.mate1 ]; then
        gzip -c star_out/Unmapped.out.mate1 > ${sample}.${index}.genome_alignment.unaligned.fastq.gz
    else
        echo -n | gzip > ${sample}.${index}.genome_alignment.unaligned.fastq.gz
    fi

    cp star_out/Log.final.out ${sample}.${index}.genome_alignment.log
    samtools idxstats -@ ${task.cpus} ${sample}.${index}.genome_alignment.bam \\
        > ${sample}.${index}.genome_alignment.stats
    echo "sample,index,program,n_reads" > ${sample}.${index}.genome_alignment.csv
    echo "${sample},${index},star," >> ${sample}.${index}.genome_alignment.csv
    """
}

GENOME_ALIGNMENT_ALIGNED.into {
    GENOME_ALIGNMENT_ALIGNED_MERGE
    GENOME_ALIGNMENT_ALIGNED_FASTQ_FASTQC }

GENOME_ALIGNMENT_UNALIGNED.into {
    GENOME_ALIGNMENT_UNALIGNED_MERGE
    GENOME_ALIGNMENT_UNALIGNED_FASTQ_FASTQC }

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
/* MERGE GENOME ALIGNMENT */
GENOME_ALIGNMENT_LOG.into { GENOME_ALIGNMENT_LOG_MERGE; GENOME_ALIGNMENT_LOG_TABLE; GENOME_ALIGNMENT_LOG_STATS }
GENOME_ALIGNMENT_CSV.into { GENOME_ALIGNMENT_CSV_INDIVIDUAL; GENOME_ALIGNMENT_CSV_MERGE }

GENOME_ALIGNMENT_LOG_TABLE
.map { sample, index, genome_log -> [ [sample, index], genome_log ] }
.set { GENOME_ALIGNMENT_LOG_TABLE_INDEXED }

GENOME_ALIGNMENT_BAM.into {
    GENOME_ALIGNMENT_BAM_FOR_BED
    GENOME_ALIGNMENT_BAM_FOR_QPASS
    GENOME_ALIGNMENT_BAM_FOR_MERGE
}

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
        set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.bam") into GENOME_ALIGNMENT_QPASS_BAM
        set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.bam.bai") into GENOME_ALIGNMENT_QPASS_BAI
        set val(sample), val(index), file("${sample}.${index}.genome.qpass.count") into GENOME_QPASS_COUNTS_OUT
        set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.stats") into GENOME_QPASS_STATS

    """
    # -F ribo_filter_flags: drop unmapped (4) + supplementary (2048); KEEP secondary (256) so
    # multi-mapper alignments contribute to every locus in the final bigwig.
    samtools view -@ ${task.cpus} -bq ${params.mapping_quality_cutoff} -F ${params.ribo_filter_flags} ${bam} > ${sample}.${index}.genome_alignment.qpass.bam && \
    samtools index -@ ${task.cpus} ${sample}.${index}.genome_alignment.qpass.bam
    samtools view -@ ${task.cpus} -c ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome.qpass.count
    samtools flagstat -@ ${task.cpus} ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome_alignment.qpass.stats
    """
}

GENOME_QPASS_COUNTS_OUT.into {
    GENOME_QPASS_COUNTS
    GENOME_QPASS_COUNTS_FOR_STATS
    GENOME_QPASS_COUNTS_FOR_NONE_DEDUP
}

GENOME_ALIGNMENT_QPASS_BAM.into {
    GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE
    GENOME_ALIGNMENT_QPASS_BAM_FOR_STATS
}
GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE.map { sample, index, bam -> [sample, bam] }.groupTuple()
.set { GENOME_ALIGNMENT_QPASS_GROUPED_BAM }

process merge_genome_bam_post_qpass {
    storeDir get_storedir('genome_alignment') + '/' + params.output.merged_lane_directory

  input:
        set val(sample), file(bam_files) from GENOME_ALIGNMENT_QPASS_GROUPED_BAM

  output:
        set val(sample), file("${sample}.genome.qpass.merged.bam"),\
                   file("${sample}.genome.qpass.merged.bam.bai") \
      into GENOME_MERGED_BAM_QPASS_PRE_DEDUP_RAW
  output:
        set val(sample), file("${sample}.genome.qpass.merged.none.bam"),\
                   file("${sample}.genome.qpass.merged.none.bam.bai") \
      into GENOME_MERGED_BAM_QPASS_NONE

  when:
    dedup_method == 'umicollapse' || dedup_method == 'none' || dedup_method == 'position'

        """
  if [ "${dedup_method}" == "umicollapse" ] || [ "${dedup_method}" == "position" ]; then
    samtools merge -@ ${task.cpus} ${sample}.genome.qpass.merged.bam ${bam_files} && samtools index -@ ${task.cpus} ${sample}.genome.qpass.merged.bam

    touch ${sample}.genome.qpass.merged.none.bam && touch ${sample}.genome.qpass.merged.none.bam.bai
  else

    touch ${sample}.genome.qpass.merged.bam && touch ${sample}.genome.qpass.merged.bam.bai
    samtools merge -@ ${task.cpus} ${sample}.genome.qpass.merged.none.bam ${bam_files} && samtools index -@ ${task.cpus} ${sample}.genome.qpass.merged.none.bam
  fi
  """
}

GENOME_MERGED_BAM_QPASS_PRE_DEDUP_RAW.set { GENOME_MERGED_BAM_QPASS_PRE_DEDUP_MAIN }

GENOME_MERGED_BAM_QPASS_PRE_DEDUP_MAIN.into {
    GENOME_MERGED_BAM_FOR_UMICOLLAPSE_DEDUP
    GENOME_MERGED_BAM_FOR_POSITION_BAM_CONV
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

    """
    samtools sort -@ ${task.cpus} -o sorted.bam ${bam}
    samtools view -@ ${task.cpus} -bq ${params.mapping_quality_cutoff} -F ${params.ribo_filter_flags} sorted.bam \
        > ${sample}.${index}.transcriptome_alignment.qpass.bam
    samtools index -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.qpass.bam
    rm sorted.bam
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
    python3 ${workflow.projectDir}/scripts/extract_reads_from_dedup_bed.py \\
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

process individual_genome_bam_to_bed {
    storeDir get_storedir('bam_to_bed') + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bam) from GENOME_ALIGNMENT_BAM_FOR_BED

    output:
        set val(sample), val(index), file("${sample}.${index}.genome.bed") into GENOME_INDIVIDUAL_BED_PRE_SPLIT
        set val(sample), val(index), file("${sample}.${index}.genome_nodedup_count.txt") \
       into GENOME_INDIVIDUAL_NODEDUP_COUNT

        """
    if [ `samtools view -@ ${task.cpus} -c ${bam}` -eq 0 ];
    then
       touch ${sample}.${index}.genome.bed
    else
        bamToBed -i ${bam} > ${sample}.${index}.genome.bed
    fi

    wc -l ${sample}.${index}.genome.bed > ${sample}.${index}.genome_nodedup_count.txt
    """
}
if (dedup_method == 'position') {
    GENOME_INDIVIDUAL_BED_PRE_SPLIT.into {
        GENOME_INDIVIDUAL_BED_FOR_INDEX
        GENOME_INDIVIDUAL_BED_FOR_SPLITTING
        GENOME_INDIVIDUAL_BED_FOR_POSITION_SPLITTING
    }
} else {
    GENOME_INDIVIDUAL_BED_PRE_SPLIT.into {
        GENOME_INDIVIDUAL_BED_FOR_INDEX
        GENOME_INDIVIDUAL_BED_FOR_SPLITTING
    }
}

process add_sample_index_col_to_genome_bed {
    storeDir get_storedir('bam_to_bed') + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bed) from GENOME_INDIVIDUAL_BED_FOR_INDEX

    output:
        set val(sample), file("${sample}.${index}.genome.with_sample_index.bed") \
         into GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED

        """
    awk -v this_sample=${sample}.${index} \
     '{ print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"this_sample) }' ${bed} \
          > ${sample}.${index}.genome.with_sample_index.bed
    """
}
GENOME_ALIGNMENT_BAM_FOR_MERGE.map { sample, index, bam -> [sample, bam] }.groupTuple()
.set { GENOME_ALIGNMENT_GROUPED_BAM }

GENOME_ALIGNMENT_ALIGNED_MERGE.map { sample, index, fastq -> [sample, fastq] }.groupTuple()
.set { GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ }

GENOME_ALIGNMENT_UNALIGNED_MERGE.map { sample, index, fastq -> [sample, fastq] }.groupTuple()
.set { GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ }
GENOME_ALIGNMENT_LOG_MERGE.map { sample, index, log -> [sample, log] }.groupTuple()
.set { GENOME_ALIGNMENT_GROUPED_LOG }

GENOME_ALIGNMENT_GROUPED_BAM.join(GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ)
                        .join(GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ)
                        .join(GENOME_ALIGNMENT_GROUPED_LOG)
                        .set { GENOME_ALIGNMENT_GROUPED_JOINT }

process merge_genome_alignment {
    storeDir get_storedir('genome_alignment') + '/' + params.output.merged_lane_directory

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

    input:

        set val(sample), file(bam), file(aligned_fastq), \
          file(unaligned_fastq), file(alignment_log) from GENOME_ALIGNMENT_GROUPED_JOINT

    output:
        set val(sample), file("${sample}.genome.bam") \
                      into GENOME_ALIGNMENT_MERGED_BAM
        set val(sample), file("${sample}.genome.bam.bai") \
                      into GENOME_ALIGNMENT_MERGED_BAI
        set val(sample), file("${sample}.genome.aligned.fastq.gz") \
                      into GENOME_ALIGNMENT_MERGED_ALIGNED_FASTQ
        set val(sample), file("${sample}.genome.unaligned.fastq.gz") \
                      into GENOME_ALIGNMENT_MERGED_UNALIGNED_FASTQ
        set val(sample), file("${sample}.genome.log") \
                      into GENOME_ALIGNMENT_MERGED_LOG
        set val(sample), file("${sample}.genome.csv") \
                      into GENOME_ALIGNMENT_MERGED_CSV

        """
    samtools merge -@ ${task.cpus} ${sample}.genome.bam ${bam} && samtools index -@ ${task.cpus} ${sample}.genome.bam && \
    (for f in ${aligned_fastq}; do if [ -s "\$f" ] && gzip -t "\$f" 2>/dev/null; then zcat "\$f"; fi; done) | gzip -c > ${sample}.genome.aligned.fastq.gz && \
    (for f in ${unaligned_fastq}; do if [ -s "\$f" ] && gzip -t "\$f" 2>/dev/null; then zcat "\$f"; fi; done) | gzip -c > ${sample}.genome.unaligned.fastq.gz && \
    cat ${alignment_log} > ${sample}.genome.log && \
    echo "sample,index,program,n_reads" > ${sample}.genome.csv && \
    echo "${sample},merged,star," >> ${sample}.genome.csv
    """
}

///////////////////////////////////////////////////////////////////////////////
//          G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////

GENOME_ALIGNMENT_MERGED_BAM.set { GENOME_BAM_FOR_DEDUP }

process genome_bam_to_bed {
    storeDir get_storedir('bam_to_bed') + '/' + params.output.individual_lane_directory

    input:
        set val(sample), file(bam) from GENOME_BAM_FOR_DEDUP

    output:
        set val(sample), file("${sample}.genome.bed") into GENOME_BED_FOR_DEDUP_PRE
        set val(sample), file("${sample}.genome_nodedup_count.txt") \
       into GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP

    when:
    dedup_method != 'none'

        """
    if [ `samtools view -@ ${task.cpus} -c ${bam}` -eq 0 ];
    then
       touch ${sample}.genome.bed
    else
        bamToBed -i ${bam} > ${sample}.genome.bed
    fi

    wc -l ${sample}.genome.bed > ${sample}.genome_nodedup_count.txt
    """
}

GENOME_BED_FOR_DEDUP_PRE.into {
    GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION
    GENOME_BED_FOR_DEDUP_PRE_FOR_STATS
}

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
                into GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION_WITH_INDEX

        when:
        dedup_method == 'position'

        """
        cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.genome.merged.pre_dedup.bed
        """
    }

    GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION_WITH_INDEX.set { GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION }
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

Channel.empty().set { GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS }

if (dedup_method == 'position') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_SEPARATION
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_BAM
    }

    process generate_post_dedup_counts_from_bed_position {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(post_dedup_bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS

        output:
        set val(sample), file("${sample}.count_after_dedup.txt") into GENOME_POSITION_DEDUP_COUNTS

        when:
        dedup_method == 'position'

        """
        wc -l < ${post_dedup_bed} > ${sample}.count_after_dedup.txt
        """
    }

    GENOME_POSITION_DEDUP_COUNTS.into {
        GENOME_POSITION_DEDUP_COUNTS_FOR_SPLITTING_RAW
        GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS_TEMP
    }

    GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS
        .mix(GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS_TEMP)
        .set { GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS }

    GENOME_POSITION_DEDUP_COUNTS_FOR_SPLITTING_RAW
        .map { sample, count_file -> [sample, count_file] }
        .set { GENOME_POSITION_DEDUP_COUNTS_BY_SAMPLE }

    GENOME_QPASS_COUNTS_FOR_STATS.into {
        GENOME_QPASS_COUNTS_FOR_STATS_SPLITTING
        GENOME_QPASS_COUNTS_FOR_STATS_SEPARATION
    }

    GENOME_POSITION_DEDUP_COUNTS_BY_SAMPLE
        .join(GENOME_QPASS_COUNTS_FOR_STATS_SPLITTING.map { sample, index, count_file -> [sample, index] }.unique(), by: 0)
        .map { sample, count_file, index -> [sample, index, count_file] }
        .set { GENOME_POSITION_DEDUP_COUNTS_FOR_SPLITTING }

    process split_post_dedup_counts_for_position {
        input:
        set val(sample), val(index), file(count_file) from GENOME_POSITION_DEDUP_COUNTS_FOR_SPLITTING

        output:
        set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING_POSITION

        when:
        dedup_method == 'position'

        """
        # For position dedup, we use the same merged count for all indices since dedup is done after merging
        # Get the count from the merged file and create count file for this specific index
        count=\$(cat ${count_file})
        echo "\${count}" > ${sample}.${index}.count_after_dedup.txt
        """
    }

    GENOME_QPASS_COUNTS_FOR_STATS_SEPARATION
        .map { sample, index, count_file -> [sample, index] }
        .unique()
        .combine(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_SEPARATION, by: 0)
        .set { GENOME_BED_FOR_POSITION_SEPARATION }

    process separate_genome_bed_post_dedup {
        storeDir get_publishdir('alignments') + '/ribo/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(bed) from GENOME_BED_FOR_POSITION_SEPARATION

        output:
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") \
                into GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION

        when:
        dedup_method == 'position'

        """
        # Extract reads for this specific lane using sample index column (column 7)
        awk -v this_sample=${sample}.${index} '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} > ${sample}.${index}.genome.post_dedup.bed
        """
    }

    process count_individual_separated_bed_position {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(bed) from GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION

        output:
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION

        when:
        dedup_method == 'position'

        """
        wc -l < ${bed} > ${sample}.${index}.count_after_dedup.txt
        """
    }
} else if (dedup_method == 'umicollapse') {
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX }
} else {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX }
}

process genome_deduplicate_umicollapse {
    storeDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam), file(bai) from GENOME_MERGED_BAM_FOR_UMICOLLAPSE_DEDUP

    output:
    set val(sample), file("${sample}.dedup.bam") into GENOME_UMI_TOOLS_DEDUP_BAM
    set val(sample), file("${sample}.dedup.bam") into GENOME_UMI_DEDUP_BAM_FOR_BED

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
    """
}

process genome_umi_dedup_bam_to_bed {
    storeDir get_publishdir('alignments') + '/ribo/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam) from GENOME_UMI_DEDUP_BAM_FOR_BED

    output:
    set val(sample), file("${sample}.genome.post_dedup.bed") into GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI

    when:
    dedup_method == 'umicollapse'

    """
    bamToBed -i ${bam} > ${sample}.genome.post_dedup.bed
    """
}

if (dedup_method == 'umicollapse') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_INDEX
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_STATS
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_MIX
    }
} else {
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_MIX }
}

GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX.mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_MIX)
.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP }

if (dedup_method == 'umicollapse') {
    GENOME_UMI_TOOLS_DEDUP_BAM.into {
        GENOME_UMI_TOOLS_DEDUP_BAM_FOR_SPLITTING
        GENOME_UMI_TOOLS_DEDUP_BAM_FOR_OTHER
    }
}

///////////////////////////////////////////////////////////////////////////////
//          P - S I T E   O F F S E T   C O R R E C T I O N
///////////////////////////////////////////////////////////////////////////////

if (dedup_method == 'umicollapse') {
    GENOME_UMI_TOOLS_DEDUP_BAM_FOR_SPLITTING
        .cross(GENOME_INDIVIDUAL_BED_FOR_SPLITTING.map { sample, index, bed -> [sample, index] }.unique())
        .map { merged_data, individual_data ->
            [individual_data[0], individual_data[1], merged_data[1]]
        }
        .set { GENOME_DEDUP_BAM_FOR_SPLITTING }

    process split_genome_dedup_bam_to_individual {
        storeDir get_publishdir('alignments') + '/ribo/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(merged_bam) from GENOME_DEDUP_BAM_FOR_SPLITTING

        output:
            set val(sample), val(index), file("${sample}.${index}.post_dedup.bam"), \
                             file("${sample}.${index}.post_dedup.bam.bai") into GENOME_INDIVIDUAL_DEDUP_BAM
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING_UMI
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") into GENOME_INDIVIDUAL_POST_DEDUP_BED_UMI

        when:
        dedup_method == 'umicollapse'

        """
        # Check if read groups exist in the merged BAM
        if ! samtools view -H ${merged_bam} | grep -q "^@RG"; then
            echo "ERROR: No read groups found in merged BAM ${merged_bam}" >&2
            echo "Read groups are required for splitting merged UMI dedup BAM" >&2
            exit 1
        fi

        # Split by read group ID
        samtools view -@ ${task.cpus} -B -r ${sample}.${index} ${merged_bam} -o ${sample}.${index}.post_dedup.bam

        # Check if splitting produced any reads
        read_count=\$(samtools view -@ ${task.cpus} -c ${sample}.${index}.post_dedup.bam)

        if [ \${read_count} -eq 0 ]; then
            echo "WARNING: No reads found for read group ${sample}.${index}" >&2
        fi

        echo \${read_count} > ${sample}.${index}.count_after_dedup.txt
        samtools index -@ ${task.cpus} ${sample}.${index}.post_dedup.bam

        if [ \${read_count} -eq 0 ]; then
            touch ${sample}.${index}.genome.post_dedup.bed
        else
            bamToBed -i ${sample}.${index}.post_dedup.bam > ${sample}.${index}.genome.post_dedup.bed
        fi
        """
    }

    GENOME_INDIVIDUAL_DEDUP_BAM.set { GENOME_INDIVIDUAL_DEDUP_BAM_FOR_STATS }
} else {
    Channel.empty().set { GENOME_INDIVIDUAL_DEDUP_BAM_FOR_STATS }
}

if (dedup_method == 'position') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_BAM
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

        when:
        dedup_method == 'position'

        """
        set +u
        source \$(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome
        set -u
        python3 ${workflow.projectDir}/scripts/extract_reads_from_dedup_bed.py \
            --bam ${qpass_bam} \
            --bed ${dedup_bed} \
            --output ${sample}.post_dedup.bam && \
        samtools index -@ ${task.cpus} ${sample}.post_dedup.bam
        """
    }

    GENOME_POST_DEDUP_BAM_POSITION_MERGED.set { GENOME_BAM_FOR_BIGWIG_FINAL }
} else if (dedup_method == 'umicollapse') {
    GENOME_UMI_TOOLS_DEDUP_BAM_FOR_OTHER.into {
        GENOME_BAM_FOR_DOWNSTREAM
        GENOME_DEDUP_BAM_FOR_INDEXING
    }

    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT }

    process generate_merged_dedup_count_umi {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(post_dedup_bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT

        output:
            set val(sample), file("${sample}.count_after_dedup.txt") into GENOME_MERGED_DEDUP_COUNT_UMI

        when:
        dedup_method == 'umicollapse'

        """
        wc -l < ${post_dedup_bed} > ${sample}.count_after_dedup.txt
        """
    }

    process index_dedup_bam_for_bigwig {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(bam) from GENOME_DEDUP_BAM_FOR_INDEXING

        output:
            set val(sample), file(bam), file("${bam}.bai") into GENOME_BAM_FOR_BIGWIG_FINAL

        when:
        dedup_method == 'umicollapse'

            """
        samtools index -@ ${task.cpus} ${bam}
        """
    }
} else {
    GENOME_MERGED_BAM_QPASS_NONE.set { GENOME_BAM_FOR_BIGWIG_FINAL }
}

if (dedup_method == 'umicollapse') {
    GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING_UMI.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
}
else if (dedup_method == 'position') {
    GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
}
else {
    GENOME_QPASS_COUNTS_FOR_NONE_DEDUP.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
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
    .set { GENOME_BAM_FOR_BIGWIG_WITH_STRAND }

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
    """
    if [ "${strand_arg}" == "F" ] || [ "${strand_arg}" == "FR" ]; then
        bamCoverage -b ${bam} -o ${sample}.ribo.plus.bigWig \
            --filterRNAstrand reverse --binSize 1 -p ${task.cpus} --minMappingQuality 0 --outFileFormat bigwig
        bamCoverage -b ${bam} -o ${sample}.ribo.minus.bigWig \
            --filterRNAstrand forward --binSize 1 -p ${task.cpus} --minMappingQuality 0 --outFileFormat bigwig
    else
        bamCoverage -b ${bam} -o ${sample}.ribo.plus.bigWig \
            --filterRNAstrand forward --binSize 1 -p ${task.cpus} --minMappingQuality 0 --outFileFormat bigwig
        bamCoverage -b ${bam} -o ${sample}.ribo.minus.bigWig \
            --filterRNAstrand reverse --binSize 1 -p ${task.cpus} --minMappingQuality 0 --outFileFormat bigwig
    fi
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

GENOME_QPASS_COUNTS
  .map { sample, index, count_file -> [ [sample, index], count_file] }
  .set { GENOME_QPASS_COUNTS_INDEXED }

GENOME_INDIVIDUAL_DEDUP_COUNT
  .map { sample, index, count_file -> [ [sample, index], count_file] }
  .set { GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED }

CLIP_LOG_INDEXED_FOR_GENOME.join(FILTER_LOG_INDEXED_FOR_GENOME)
            .join(GENOME_ALIGNMENT_LOG_STATS_INDEXED)
            .join(GENOME_QPASS_COUNTS_INDEXED)
            .join(GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED)
            .map { key, clip_log, filter_log, genome_log, qpass_count, dedup_count ->
                [key[0], key[1], clip_log, filter_log, genome_log, qpass_count, dedup_count]
            }
            .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }

process individual_genome_alignment_stats {
    storeDir get_storedir('stats') + '/genome/individual'

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

input:
    set val(sample), val(index), file(clip_log), file(filter_log),\
    file(genome_log), file(qpass_count),\
    file(dedup_count)\
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

# STAR genome alignment log — parsed by scripts/parse_star_log.py
import csv, subprocess
proc = subprocess.run(
    ['python3', '${baseDir}/scripts/parse_star_log.py', '${genome_log}'],
    check=True, capture_output=True, text=True,
)
star_row = next(csv.DictReader(proc.stdout.splitlines(), delimiter='\\t'))
genome_once  = int(star_row['uniquely_mapped'])
genome_multi = int(star_row['multi_loci_mapped'])
genome_unal  = int(star_row['unmapped_total'])
genome_total = int(star_row['primary_aligned_total'])

with open('${qpass_count}') as fh:
    qpass = int(fh.read().strip().split()[0])
with open('${dedup_count}') as fh:
    after_dedup = int(fh.read().strip().split()[0])

rows = [('total_reads', total_reads), ('clipped_reads', clipped_reads),
        ('filtered_out', filtered_out), ('filter_kept', filter_kept),
        ('genome_aligned_once', genome_once), ('genome_aligned_many', genome_multi),
        ('genome_total_aligned', genome_total), ('genome_unaligned', genome_unal),
        ('genome_qpass_aligned_reads', qpass), ('genome_after_dedup', after_dedup)]
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
        echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> genome_individual_essential.csv
        '''
    } else {
        """
  rfc merge overall-stats \
   -o raw_combined_individual_genome_aln_stats.csv \
      ${stat_table} && \
  rfc stats-percentage \
  -i raw_combined_individual_genome_aln_stats.csv \
  -l genome \
  -o genome_individual_essential.csv
        """
    }
}

GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GROUPING
  .map { sample, index, file -> [ sample, file ] }
  .groupTuple()
  .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED }

if (dedup_method == 'position') {
    GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
        .combine(GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS, by: 0)
        .map { sample, stat_files, merged_dedup_count_file ->
            [sample, stat_files, merged_dedup_count_file]
        }
        .set { GENOME_STATS_INPUT }
} else if (dedup_method == 'umicollapse') {
    GENOME_MERGED_DEDUP_COUNT_UMI
        .map { sample, count_file -> [sample, count_file] }
        .set { GENOME_UMI_DEDUP_COUNTS_FOR_MERGED_STATS }
    GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
        .combine(GENOME_UMI_DEDUP_COUNTS_FOR_MERGED_STATS, by: 0)
        .map { sample, stat_files, merged_dedup_count_file ->
            [sample, stat_files, merged_dedup_count_file]
        }
        .set { GENOME_STATS_INPUT }
} else {
    GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
        .map { sample, stat_files ->
            [sample, stat_files, null]
        }
        .set { GENOME_STATS_INPUT }
}

process sum_individual_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/merged'

    beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

input:
    set val(sample), file(stat_files), val(merged_dedup_count_file) from GENOME_STATS_INPUT

output:
    set val(sample), file("${sample}.genome_merged.csv") into GENOME_MERGED_ALIGNMENT_STATS

    script:
    merged_count_file_exists = merged_dedup_count_file != null && merged_dedup_count_file.toString() != 'null' && merged_dedup_count_file.toString() != ''

    if ((dedup_method == 'position' || dedup_method == 'umicollapse') && merged_count_file_exists) {
        """
        rfc sum-stats -n ${sample} -o ${sample}.genome_merged.tmp.csv ${stat_files}

        python3 ${workflow.projectDir}/scripts/update_merged_stats_with_counts.py \\
            --dedup-count-file ${merged_dedup_count_file} \\
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
        // Create empty stats file when no input is available
        '''
        echo "No merged genome statistics data available" > genome_merged_essential.csv
        echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> genome_merged_essential.csv
        '''
    } else {
        """
    rfc merge overall-stats \
    -o raw_combined_merged_genome_aln_stats.csv \
    ${stat_files} && \
    rfc stats-percentage \
    -i raw_combined_merged_genome_aln_stats.csv \
    -l genome \
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

    // RNA-seq is paired-end only. YAML schema:
    //   [R1, R2] list -> paired-end
    // Channel tuple shape: [sample, index, file(r1), file(r2)].
    def rnaseq_lane_tuples = []
    params.get('rnaseq', [:]).fastq.each { sample_name, lanes ->
        lanes.eachWithIndex { entry, i ->
            def lane_idx = i + 1
            if ((entry instanceof List) && entry.size() == 2) {
                rnaseq_lane_tuples << [sample_name, lane_idx,
                    file("${rnaseq_fastq_base}${entry[0]}"),
                    file("${rnaseq_fastq_base}${entry[1]}")]
            } else {
                throw new IllegalArgumentException(
                    "rnaseq.fastq.${sample_name} entry #${lane_idx} must be a paired-end list [R1, R2]. Got: ${entry}")
            }
        }
    }

    Channel.from(rnaseq_lane_tuples)
    .into {  RNASEQ_FASTQ
             RNASEQ_FASTQ_VERBOSE
             RNASEQ_FASTQ_FASTQC
             RNASEQ_FASTQ_CLIP
             RNASEQ_FASTQ_EXISTENCE}

    if (params.do_check_file_existence) {
        RNASEQ_FASTQ_EXISTENCE
            .map { sample, index, r1, r2 ->
                file_exists(r1) && file_exists(r2)
            }
    }

    process rnaseq_raw_fastqc {
        publishDir get_publishdir('fastqc', true), mode: 'copy'

      input:
        set val(sample), val(index), file(r1), file(r2) from RNASEQ_FASTQ_FASTQC

      output:
        file("*_fastqc.html") into RNASEQ_FASTQC_HTML
        file("*_fastqc.zip")  into RNASEQ_FASTQC_OUT

    when:
    params.do_fastqc && do_rnaseq

        script:
        """
        ln -sf ${r1} ${sample}.${index}.R1.fastq.gz
        ln -sf ${r2} ${sample}.${index}.R2.fastq.gz
        fastqc ${sample}.${index}.R1.fastq.gz ${sample}.${index}.R2.fastq.gz --outdir=\$PWD -t ${task.cpus}
        """
    }

    process rnaseq_clip {
        storeDir get_storedir('clip', true)

  input:
        set val(sample), val(index), file(r1), file(r2) from RNASEQ_FASTQ_CLIP

  output:
        set val(sample), val(index),
            file("${sample}.${index}.clipped_R1.fastq.gz"),
            file("${sample}.${index}.clipped_R2.fastq.gz") into RNASEQ_CLIP_OUT
        set val(sample), val(index), file("${sample}.${index}.clipped.log") \
                                                      into RNASEQ_CLIP_LOG

        script:
        """
        cutadapt --cores=${task.cpus} ${params.rnaseq.clip_arguments} \\
            -o ${sample}.${index}.clipped_R1.fastq.gz \\
            -p ${sample}.${index}.clipped_R2.fastq.gz \\
            ${r1} ${r2} \\
            > ${sample}.${index}.clipped.log 2>&1
        """
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
        set val(sample), val(index), file(r1), file(r2) from RNASEQ_CLIP_OUT
        set val(bowtie2_index_base), file(bowtie2_index_files) \
                          from RNASEQ_FILTER_INDEX.first()

    output:
        set val(sample), val(index), file("${sample}.${index}.filter.bam") \
        into RNASEQ_FILTER_BAM
        set val(sample), val(index), file("${sample}.${index}.filter.bam.bai") \
        into RNASEQ_FILTER_BAI
        set val(sample), val(index),
            file("${sample}.${index}.aligned.filter_R1.fastq.gz"),
            file("${sample}.${index}.aligned.filter_R2.fastq.gz") into RNASEQ_FILTER_ALIGNED
        set val(sample), val(index),
            file("${sample}.${index}.unaligned.filter_R1.fastq.gz"),
            file("${sample}.${index}.unaligned.filter_R2.fastq.gz") into RNASEQ_FILTER_UNALIGNED
        set val(sample), val(index), file("${sample}.${index}.filter.log") \
        into RNASEQ_FILTER_LOG
        set val(sample), val(index), file("${sample}.${index}.filter.stats") \
        into RNASEQ_FILTER_STATS

        // bowtie2 --un-conc-gz/--al-conc-gz expand `%` to 1/2, so passing
        // <prefix>_R%.fastq.gz produces <prefix>_R1.fastq.gz and <prefix>_R2.fastq.gz
        // directly. A pair lands in "unaligned" only when it does NOT align
        // concordantly to the rRNA index.
        script:
        """
        set -o pipefail
        bowtie2 ${params.rnaseq.filter_arguments} \\
                -x ${bowtie2_index_base} -1 ${r1} -2 ${r2} \\
                --threads ${task.cpus} \\
                --al-conc-gz ${sample}.${index}.aligned.filter_R%.fastq.gz \\
                --un-conc-gz ${sample}.${index}.unaligned.filter_R%.fastq.gz \\
                2> ${sample}.${index}.filter.log \\
                | samtools view -@ ${task.cpus} -b - \\
                | samtools sort -@ ${task.cpus} -o ${sample}.${index}.filter.bam \\
                && samtools index -@ ${task.cpus} ${sample}.${index}.filter.bam \\
                && samtools idxstats -@ ${task.cpus} ${sample}.${index}.filter.bam  > \\
                   ${sample}.${index}.filter.stats
        """
    }

    // Set RNA-seq filter log for genome stats
    RNASEQ_FILTER_LOG.set { RNASEQ_FILTER_LOG_FOR_GENOME }

    RNASEQ_FILTER_UNALIGNED.into { RNASEQ_FILTER_UNALIGNED_FASTQ_READ_LENGTH
        RNASEQ_FILTER_UNALIGNED_FASTQ_FASTQC
        RNASEQ_FILTER_UNALIGNED_GENOME}

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
        set val(sample), val(index), file(r1), file(r2) from RNASEQ_GENOME_INPUT_CHANNEL
        set val(genome_base), file(genome_files) from RNASEQ_GENOME_INDEX_FOR_ALIGN.first()

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.bam") \
        into RNASEQ_GENOME_ALIGNMENT_BAM
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.bam.bai") \
        into RNASEQ_GENOME_ALIGNMENT_BAI
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.aligned.fastq.gz") \
        into RNASEQ_GENOME_ALIGNMENT_ALIGNED
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.unaligned.fastq.gz") \
        into RNASEQ_GENOME_ALIGNMENT_UNALIGNED
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.log") \
        into RNASEQ_GENOME_ALIGNMENT_LOG
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.stats") \
        into RNASEQ_GENOME_ALIGNMENT_STATS

        // Mirrors ENCODE rna-seq-pipeline alignment (align.py).
        // --outSAMunmapped Within keeps unmapped reads in the BAM. Paired-end:
        // R1 and R2 are passed to STAR as space-separated readFiles.
        script:
        """
    set -o pipefail
    mkdir -p star_out

    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${r1} ${r2} \\
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

    # Both fastq slots are unused downstream — emit empty files for channel compatibility.
    echo -n | gzip > ${sample}.${index}.rnaseq_genome_alignment.aligned.fastq.gz
    echo -n | gzip > ${sample}.${index}.rnaseq_genome_alignment.unaligned.fastq.gz

    cp star_out/Log.final.out ${sample}.${index}.rnaseq_genome_alignment.log
    samtools idxstats -@ ${task.cpus} ${sample}.${index}.rnaseq_genome_alignment.bam \\
        > ${sample}.${index}.rnaseq_genome_alignment.stats
    echo "sample,index,program,n_reads" > ${sample}.${index}.rnaseq_genome_alignment.csv
    echo "${sample},${index},star," >> ${sample}.${index}.rnaseq_genome_alignment.csv
        """
    }

    // Split RNA-seq genome channels for different uses
    RNASEQ_GENOME_ALIGNMENT_BAM.into {
        RNASEQ_GENOME_ALIGNMENT_BAM_FOR_QUALITY
        RNASEQ_GENOME_ALIGNMENT_BAM_FOR_MERGE
    }

    RNASEQ_GENOME_ALIGNMENT_LOG.set { RNASEQ_GENOME_ALIGNMENT_LOG_FOR_STATS }

    // RNA-seq genome quality filtering
    process rnaseq_genome_quality_filter {
        storeDir get_storedir('quality_filter', true)

    input:
        set val(sample), val(index), file(bam) from RNASEQ_GENOME_ALIGNMENT_BAM_FOR_QUALITY

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.qpass.bam") into RNASEQ_GENOME_ALIGNMENT_QPASS_BAM
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.qpass.bam.bai") into RNASEQ_GENOME_ALIGNMENT_QPASS_BAI
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.qpass.count") into RNASEQ_GENOME_QPASS_COUNTS
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.qpass.stats") into RNASEQ_GENOME_QPASS_STATS

        script:
        def rnaseq_cfg   = params.get('rnaseq', [:])
        def rna_mapq     = rnaseq_cfg.get('mapping_quality_cutoff', 4)
        def rna_flags    = rnaseq_cfg.get('filter_flags', 2308)
        """
    # -F rnaseq.filter_flags: drop unmapped (4) + secondary (256) + supplementary (2048).
    samtools view -@ ${task.cpus} -bq ${rna_mapq} -F ${rna_flags} ${bam} \
        > ${sample}.${index}.rnaseq_genome_alignment.qpass.bam && \\
    samtools index -@ ${task.cpus} ${sample}.${index}.rnaseq_genome_alignment.qpass.bam
    samtools view -@ ${task.cpus} -c ${sample}.${index}.rnaseq_genome_alignment.qpass.bam \
        > ${sample}.${index}.rnaseq_genome.qpass.count
    samtools flagstat -@ ${task.cpus} ${sample}.${index}.rnaseq_genome_alignment.qpass.bam \
        > ${sample}.${index}.rnaseq_genome_alignment.qpass.stats
    """
    }

    // Split quality-filtered BAM channel for BED conversion and nodedup bigwig generation
    RNASEQ_GENOME_ALIGNMENT_QPASS_BAM.into {
        RNASEQ_GENOME_QPASS_BAM_FOR_BED
        RNASEQ_GENOME_QPASS_BAM_FOR_NODEDUP_MERGE
    }

    // Convert individual RNA-seq genome BAMs to BED
    process rnaseq_individual_genome_bam_to_bed {
        storeDir get_storedir('bam_to_bed', true) + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bam) from RNASEQ_GENOME_QPASS_BAM_FOR_BED

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.bed") into RNASEQ_INDIVIDUAL_GENOME_BED
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_nodedup_count.txt") into RNASEQ_INDIVIDUAL_GENOME_NODEDUP_COUNT

        """
    if [ `samtools view -@ ${task.cpus} -c ${bam}` -eq 0 ];
    then
       touch ${sample}.${index}.rnaseq_genome.bed
    else
        bamToBed -i ${bam} > ${sample}.${index}.rnaseq_genome.bed
    fi

    wc -l ${sample}.${index}.rnaseq_genome.bed > ${sample}.${index}.rnaseq_genome_nodedup_count.txt
    """
    }

    // Split RNA-seq genome BED channel for different uses
    RNASEQ_INDIVIDUAL_GENOME_BED.into {
        RNASEQ_GENOME_BED_FOR_INDEX_COL
        RNASEQ_GENOME_BED_FOR_STATS
        RNASEQ_GENOME_BED_FOR_INDEX_SEP_PRE
    }

    // Add sample index column to RNA-seq genome bed files
    process rnaseq_add_sample_index_col_to_genome_bed {
        storeDir get_storedir('bam_to_bed', true) + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bed) from RNASEQ_GENOME_BED_FOR_INDEX_COL

    output:
        set val(sample), file("${sample}.${index}.rnaseq_genome.with_sample_index.bed") into RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED

        """
    awk -v newcol=${sample}.${index} '{print(\$0"\\t"newcol)}' ${bed} > ${sample}.${index}.rnaseq_genome.with_sample_index.bed
    """
    }

    // Group genome BED files with sample index for merging
    RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED.groupTuple()
    .set { RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED }

    // Merge RNA-seq genome BED files)
    process rnaseq_merge_genome_bed {
        storeDir get_storedir('bam_to_bed', true) + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(bed_files) from RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED

    output:
        set val(sample), file("${sample}.merged.pre_dedup.bed") \
     into RNASEQ_GENOME_BED_MERGED_PRE_DEDUP

        """
    cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.merged.pre_dedup.bed
    """
    }

    RNASEQ_GENOME_BED_MERGED_PRE_DEDUP.set { RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP }

    // Create genome nodedup count files from individual BED files
    process rnaseq_genome_create_nodedup_counts {
    input:
        set val(sample), val(index), file(bed) from RNASEQ_GENOME_BED_FOR_STATS

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_nodedup_count.txt") into RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP

        """
    wc -l ${bed} > ${sample}.${index}.rnaseq_genome_nodedup_count.txt
    """
    }

    // RNA-seq genome individual statistics compilation
    // Create indexed channels directly and use them without splitting to avoid consumption conflicts
    RNASEQ_CLIP_LOG_FOR_GENOME.map { sample, index, clip_log -> [ [sample, index], clip_log ] }
        .set { RNASEQ_CLIP_LOG_INDEXED_FOR_GENOME }
    RNASEQ_FILTER_LOG_FOR_GENOME.map { sample, index, filter_log -> [ [sample, index], filter_log ] }
        .set { RNASEQ_FILTER_LOG_INDEXED_FOR_GENOME }
    RNASEQ_GENOME_ALIGNMENT_LOG_FOR_STATS.map { sample, index, genome_log -> [ [sample, index], genome_log ] }
        .set { RNASEQ_GENOME_ALIGNMENT_LOG_INDEXED }
    RNASEQ_GENOME_QPASS_COUNTS.map { sample, index, qpass_count -> [ [sample, index], qpass_count ] }
        .set { RNASEQ_GENOME_QPASS_COUNTS_INDEXED }

    // Group per-lane BAMs by sample for merging. The merged aligned/unaligned
    // fastqs, concatenated log, and csv that this process used to emit have no
    // downstream consumers (only the merged BAM+BAI are used), so we drop them
    // entirely — they were single-threaded gzip reencodes of gigabytes of dead data.
    RNASEQ_GENOME_ALIGNMENT_BAM_FOR_MERGE.map { sample, index, bam -> [sample, bam] }.groupTuple()
    .set { RNASEQ_GENOME_ALIGNMENT_GROUPED_BAM }

    process rnaseq_merge_genome_alignment {
        storeDir get_storedir('genome_alignment', true) + '/' + params.output.merged_lane_directory

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

    input:
        set val(sample), file(bam) from RNASEQ_GENOME_ALIGNMENT_GROUPED_BAM

    output:
        set val(sample), file("${sample}.rnaseq_genome.bam") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_BAM
        set val(sample), file("${sample}.rnaseq_genome.bam.bai") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_BAI

        """
    samtools merge -@ ${task.cpus} ${sample}.rnaseq_genome.bam ${bam}
    samtools index -@ ${task.cpus} ${sample}.rnaseq_genome.bam
    """
    }

///////////////////////////////////////////////////////////////////////////////
//          R N A - S E Q   G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////

    // Merge quality-filtered BAMs for nodedup bigwig generation
    RNASEQ_GENOME_QPASS_BAM_FOR_NODEDUP_MERGE
        .map { sample, index, bam -> [sample, bam] }
        .groupTuple()
        .set { RNASEQ_GENOME_QPASS_GROUPED_BAM }

    // Merge quality-filtered BAMs for nodedup bigwig
    process rnaseq_merge_quality_filtered_bam {
        storeDir get_storedir('genome_alignment', true) + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(bam_files) from RNASEQ_GENOME_QPASS_GROUPED_BAM

        output:
        set val(sample), file("${sample}.rnaseq_genome.qpass.merged.bam") \
            into RNASEQ_GENOME_QPASS_MERGED_BAM
        set val(sample), file("${sample}.rnaseq_genome.qpass.merged.bam.bai") \
            into RNASEQ_GENOME_QPASS_MERGED_BAI

        """
        samtools merge -@ ${task.cpus} ${sample}.rnaseq_genome.qpass.merged.bam ${bam_files} && \
        samtools index -@ ${task.cpus} ${sample}.rnaseq_genome.qpass.merged.bam
        """
    }

    RNASEQ_GENOME_QPASS_MERGED_BAM.into {
        RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_NODEDUP
        RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_DEDUP
    }

    RNASEQ_GENOME_ALIGNMENT_MERGED_BAM.set { RNASEQ_GENOME_BAM_FOR_BIGWIG_REF }
    RNASEQ_GENOME_ALIGNMENT_MERGED_BAI.set { RNASEQ_GENOME_BAI_FOR_BIGWIG_REF }

///////////////////////////////////////////////////////////////////////////////
//  R N A - S E Q   B I G W I G   G E N E R A T I O N   (N O D E D U P)
///////////////////////////////////////////////////////////////////////////////

    // Use quality-filtered merged BAM for nodedup bigwig generation
    RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_NODEDUP
        .join(RNASEQ_GENOME_QPASS_MERGED_BAI)
        .set { RNASEQ_GENOME_FOR_NODEDUP_BIGWIG }

    // RNA-seq genome deduplication now follows transcriptome pattern BED-based

    // RNA-seq genome deduplication (position-based only)
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

    // Split dedup BED channel for separation, count generation, and bigwig generation
    RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.into {
        RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_SEP
        RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_COUNT
        RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_BIGWIG
    }

    // Individual count files will be created by the separation process

    // Create channel for genome separation (following transcriptome pattern)
    RNASEQ_GENOME_BED_FOR_INDEX_SEP_PRE
    .map { sample, index, file -> [sample, index] }
    .combine(RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_SEP, by:0)
    .set { RNASEQ_GENOME_BED_FOR_INDEX_SEP_POST_DEDUP }

    // Separate RNA-seq genome merged post_dedup bed back into individual lane files
    process rnaseq_separate_genome_bed_post_dedup {
        storeDir  get_publishdir('alignments') + '/rnaseq/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bed) from RNASEQ_GENOME_BED_FOR_INDEX_SEP_POST_DEDUP

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.post_dedup.bed") \
     into RNASEQ_GENOME_BED_DEDUPLICATED
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.count_after_dedup.txt")\
     into RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP_SEPARATED_RAW

    when:
    rnaseq_dedup_method == 'position'

        """
    awk -v this_sample=${sample}.${index} \
     '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} > ${sample}.${index}.rnaseq_genome.post_dedup.bed \
      && wc -l ${sample}.${index}.rnaseq_genome.post_dedup.bed > ${sample}.${index}.rnaseq_genome.count_after_dedup.txt
    """
    }

    if (rnaseq_dedup_method == 'position') {
        RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP_SEPARATED_RAW
            .set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT }
    } else {
        RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP
            .set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT }
    }

    RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT
        .map { sample, index, dedup_count -> [ [sample, index], dedup_count ] }
        .set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED }

    RNASEQ_CLIP_LOG_INDEXED_FOR_GENOME.join(RNASEQ_FILTER_LOG_INDEXED_FOR_GENOME)
            .join(RNASEQ_GENOME_ALIGNMENT_LOG_INDEXED)
            .join(RNASEQ_GENOME_QPASS_COUNTS_INDEXED)
            .join(RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED)
            .flatten()
            .collate(7)
            .set { RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }

    process rnaseq_individual_genome_alignment_stats {
        storeDir get_storedir('stats', true) + '/genome/individual'

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

        input:
        set val(sample), val(index), file(clip_log), file(filter_log), \
            file(genome_log), file(qpass_count), file(dedup_count) \
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

# filter log (bowtie2)
fl = [l for l in open('${filter_log}') if l.strip() and not l[0].isalpha()]
filtered_out = int(fl[3].split()[0]) + int(fl[4].split()[0])
filter_kept  = int(fl[2].split()[0])

# STAR genome alignment log — parsed by scripts/parse_star_log.py
import csv, subprocess
proc = subprocess.run(
    ['python3', '${baseDir}/scripts/parse_star_log.py', '${genome_log}'],
    check=True, capture_output=True, text=True,
)
star_row = next(csv.DictReader(proc.stdout.splitlines(), delimiter='\\t'))
genome_once  = int(star_row['uniquely_mapped'])
genome_multi = int(star_row['multi_loci_mapped'])
genome_unal  = int(star_row['unmapped_total'])
genome_total = int(star_row['primary_aligned_total'])

with open('${qpass_count}') as fh:
    qpass = int(fh.read().strip().split()[0])
with open('${dedup_count}') as fh:
    after_dedup = int(fh.read().strip().split()[0])

rows = [('total_reads', total_reads), ('clipped_reads', clipped_reads),
        ('filtered_out', filtered_out), ('filter_kept', filter_kept),
        ('genome_aligned_once', genome_once), ('genome_aligned_many', genome_multi),
        ('genome_total_aligned', genome_total), ('genome_unaligned', genome_unal),
        ('genome_qpass_aligned_reads', qpass), ('genome_after_dedup', after_dedup)]
with open('${sample}.${index}.rnaseq_genome_individual.csv', 'w') as fh:
    fh.write(',${sample}.${index}\\n')
    for k, v in rows:
        fh.write(f'{k},{v}\\n')
PYEOF
        """
    }

///////////////////////////////////////////////////////////////////////////////
//  R N A - S E Q   B I G W I G   G E N E R A T I O N  (MERGED)
///////////////////////////////////////////////////////////////////////////////

    // Generate merged dedup count from merged post_dedup BED for RNA-seq position dedup stats
    process rnaseq_generate_merged_dedup_count {
        storeDir get_storedir('alignment_ribo', true) + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(post_dedup_bed) from RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_COUNT

        output:
        set val(sample), file("${sample}.rnaseq_genome.count_after_dedup.txt") \
            into RNASEQ_GENOME_MERGED_DEDUP_COUNT

        when:
        rnaseq_dedup_method == 'position'

        """
        wc -l ${post_dedup_bed} > ${sample}.rnaseq_genome.count_after_dedup.txt
        """
    }

    // Group RNA-seq dedup BED files by experiment for bigWig generation
    // Extract experiment name from sample (e.g., "experiment_1" from "experiment_1.1")
    RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_BIGWIG
        .map { sample, bed ->
            // Extract experiment name by removing the last dot and everything after it
            // This handles samples like "experiment_1.1" -> "experiment_1"
            def experiment = sample.contains('.') ? sample.substring(0, sample.lastIndexOf('.')) : sample
            [experiment, [sample, bed]]
        }
        .groupTuple(by: 0)
        .map { experiment, sample_bed_tuples ->
            // Unpack grouped tuples: each tuple is [sample, bed]
            def samples = []
            def beds = []
            sample_bed_tuples.each { tuple ->
                samples << tuple[0]
                beds << tuple[1]
            }
            [experiment, samples, beds]
        }
        .set { RNASEQ_GENOME_DEDUP_BED_GROUPED_BY_EXPERIMENT }

    // Merge BED files per experiment before bigWig generation
    process rnaseq_merge_dedup_beds_by_experiment {
        storeDir get_storedir('alignment_ribo', true) + '/' + params.output.merged_lane_directory

        input:
        set val(experiment), val(samples), file(beds) from RNASEQ_GENOME_DEDUP_BED_GROUPED_BY_EXPERIMENT

        output:
        set val(experiment), file("${experiment}.rnaseq_genome.post_dedup.merged.bed") \
            into RNASEQ_GENOME_DEDUP_BED_MERGED_FOR_BIGWIG

        when:
        rnaseq_dedup_method == 'position'

        """
        # Merge all BED files for this experiment
        cat ${beds} | sort -k1,1 -k2,2n -k3,3n > ${experiment}.rnaseq_genome.post_dedup.merged.bed
        """
    }

    RNASEQ_GENOME_DEDUP_BED_MERGED_FOR_BIGWIG
        .combine(RNASEQ_GENOME_QPASS_MERGED_BAM_FOR_DEDUP.first())
        .map { experiment, bed, sample, qpass_bam -> [experiment, bed, qpass_bam] }
        .set { RNASEQ_GENOME_DEDUP_BED_WITH_QPASS_BAM }

    process rnaseq_extract_dedup_reads_from_bam {
        storeDir get_publishdir('alignments') + '/rnaseq/' + params.output.merged_lane_directory

        input:
        set val(experiment), file(post_dedup_bed), file(qpass_bam) from RNASEQ_GENOME_DEDUP_BED_WITH_QPASS_BAM

        output:
        set val(experiment), file("${experiment}.rnaseq_genome.post_dedup.bam") \
            into RNASEQ_GENOME_DEDUP_BAM_FOR_BIGWIG
        set val(experiment), file("${experiment}.rnaseq_genome.post_dedup.bam.bai") \
            into RNASEQ_GENOME_DEDUP_BAI_FOR_BIGWIG

        when:
        rnaseq_dedup_method == 'position'

        """
        set +u
        source \$(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome
        set -u
        python3 ${workflow.projectDir}/scripts/extract_reads_from_dedup_bed.py \\
            --bam ${qpass_bam} \\
            --bed ${post_dedup_bed} \\
            --output ${experiment}.rnaseq_genome.post_dedup.bam && \\
        samtools index -@ ${task.cpus} ${experiment}.rnaseq_genome.post_dedup.bam
        """
    }

    process rnaseq_create_bigwig {
        storeDir get_publishdir('bigwigs') + '/rnaseq'

        input:
        set val(sample), file(bam), file(bai) from RNASEQ_GENOME_FOR_NODEDUP_BIGWIG

        output:
        set val(sample), file("${sample}.rnaseq.bigWig") into RNASEQ_BIGWIGS

        """
        bamCoverage -b ${bam} -o ${sample}.rnaseq.bigWig --outFileFormat bigwig \\
            --binSize 1 -p ${task.cpus}
        """
    }

///////////////////////////////////////////////////////////////////////////////
//  E N D   R N A - S E Q   B I G W I G   G E N E R A T I O N  (MERGED)
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
            echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> rnaseq_individual_stats.csv
            '''
        } else {
            """
      rfc merge overall-stats \
       -o raw_combined_individual_rnaseq_genome_aln_stats.csv \
          ${stat_table} && \
      rfc stats-percentage \
       -i raw_combined_individual_rnaseq_genome_aln_stats.csv \
       -l genome \
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
            .combine(RNASEQ_GENOME_MERGED_DEDUP_COUNT, by: 0)
            .set { RNASEQ_GENOME_STATS_INPUT }
    } else {
        RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
            .map { sample, stat_files -> [sample, stat_files, null] }
            .set { RNASEQ_GENOME_STATS_INPUT }
    }

    process sum_individual_rnaseq_genome_alignment_stats {
        executor 'local'
        storeDir get_storedir('stats', true) + '/genome/merged'

        beforeScript 'source $(conda info --base)/etc/profile.d/conda.sh && conda activate ribo_genome'

        input:
        set val(sample), file(stat_files), val(merged_dedup_count_file) from RNASEQ_GENOME_STATS_INPUT

        output:
        set val(sample), file("${sample}.rnaseq_genome_merged.csv") into RNASEQ_GENOME_MERGED_ALIGNMENT_STATS

        script:
        merged_count_file_exists = merged_dedup_count_file != null && merged_dedup_count_file.toString() != 'null' && merged_dedup_count_file.toString() != ''

        if (rnaseq_dedup_method == 'position' && merged_count_file_exists) {
            """
            rfc sum-stats -n ${sample} -o ${sample}.rnaseq_genome_merged.tmp.csv ${stat_files}

            python3 ${workflow.projectDir}/scripts/update_rnaseq_merged_stats.py \\
                --dedup-count-file ${merged_dedup_count_file} \\
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
            echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> rnaseq_stats.csv
            '''
        } else {
            """
      rfc merge overall-stats \
       -o raw_combined_merged_rnaseq_genome_aln_stats.csv \
          ${stat_table} && \
      rfc stats-percentage \
       -i raw_combined_merged_rnaseq_genome_aln_stats.csv \
       -l genome \
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
