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
  // dedup_arg: User provided dedup method
  // dedup_old: Backward compatibility boolean parameter
  // Methods: umi_tools, position, none
    def valid_methods = ['position', 'umi_tools', 'none']

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
        // Backward compatibility: default to position-based dedup
        if ( dedup_old.toLowerCase() != 'false' ) {
            return("position")
        }
    else {
            return("none")
    }
  }
}

// Command building helpers
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

// Set default input structure if not provided
if (!params.containsKey('input')) {
    params.input = [
        reference: [
            filter: './rf_sample_data/filter/human_rtRNA*',
            regions: './rf_sample_data/annotation/appris_human_24_01_2019_actual_regions.bed',
            transcript_lengths: './rf_sample_data/annotation/appris_human_24_01_2019_selected.lengths.tsv',
            genome: './rf_sample_data/genome/hisat2_spliced_index*'
        ],
        fastq: [:]
    ]
}

// Set default output structure if not provided
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
            deduplication: 'deduplication',
            bigwigs: 'bigwigs'
        ],
        output: [
            base: 'output',
            log: 'log',
            fastqc: 'fastqc',
            ribo: 'ribo'
        ],
        individual_lane_directory: 'individual',
        merged_lane_directory: 'merged'
    ]
}

// Set default boolean parameters
if (!params.containsKey('do_check_file_existence')) {
    params.do_check_file_existence = true
}
if (!params.containsKey('do_fastqc')) {
    params.do_fastqc = false
}
if (!params.containsKey('do_rnaseq')) {
    params.do_rnaseq = true
}
if (!params.containsKey('do_ribo_creation')) {
    params.do_ribo_creation = false
}
if (!params.containsKey('do_metadata')) {
    params.do_metadata = false
}

// RNA-seq deduplication method (independent from ribo-seq dedup_method)
rnaseq_dedup_method = params.get('rnaseq', [:]).get('dedup_method', 'position').toString()

fastq_base = params.get('input', [:]).get('fastq_base', '')
if (! fastq_base.endsWith('/') && fastq_base != '') {
    fastq_base = "${fastq_base}/"
}

// Group input files into [sample, fileindex, path_to_fastq_file] tuples
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

// Create correspondence log file
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

boolean hisat2_ref_exists(hisat2_ref) {
    Channel.from(['1.ht2', '2.ht2', '3.ht2', '4.ht2', '5.ht2', '6.ht2', '7.ht2', '8.ht2'])
    .map { this_suffix -> file_exists("${hisat2_ref}.${this_suffix}".replaceAll('\\*', '')) }
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
    hisat2_ref_exists(params.get('input', [:]).reference.genome)

    if (params.get('input', [:]).reference.get('post_genome', false)) {
        bt2_ref_exists(params.get('input', [:]).reference.post_genome)
    }

    file_exists(params.get('input', [:]).reference.regions)
    file_exists(params.get('input', [:]).reference.transcript_lengths)

    root_meta_file = params.get('input', [:]).get('root_meta', false)
    if (root_meta_file) {
        file_exists(root_meta_file)
    }
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
  dedup_method == 'umi_tools'

    """
  umi_tools extract -I ${fastq} -S ${sample}.${index}.umi_extracted.fastq.gz \
     -L ${sample}.${index}.umi_extracted.log \
     ${params.get('umi_tools_extract_arguments', '')}
  """
}

///////////////////////////////////////////////////////////////////////////////////////

if (dedup_method == 'none' || dedup_method == 'position') {
    CLIP_OUT_FILTER.set { FILTER_INPUT_FASTQ }
}
else if (dedup_method == 'umi_tools') {
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

// Map reads against rRNA, tRNA, adapters - keep unaligned reads for downstream

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
    bowtie2 ${params.alignment_arguments.filter} \
            -x ${bowtie2_index_base} -q ${fastq} \
            --threads ${task.cpus} \
            --al-gz ${sample}.${index}.aligned.filter.fastq.gz \
            --un-gz ${sample}.${index}.unaligned.filter.fastq.gz \
                     2> ${sample}.${index}.filter.log \
            | samtools view -bS - \
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
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* GENOME ALIGNMENT */

// Genome alignment is now the primary and required alignment method
GENOME_INPUT_CHANNEL = FILTER_UNALIGNED_GENOME

GENOME_INDEX = Channel.from([[
            params.input.reference.genome
            .split('/')[-1]
            .replaceAll('\\*$', '')
            .replaceAll('\\.$', ''),
         file(params.input.reference.genome),
        ]])

process genome_alignment {
    storeDir get_storedir('genome_alignment') + '/' + params.get('output', [:]).get('individual_lane_directory', 'individual')

input:
    set val(sample), val(index), file(fastq) from GENOME_INPUT_CHANNEL
    set val(genome_base), file(genome_files) from GENOME_INDEX.first()

output:
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.bam") \
    into GENOME_ALIGNMENT_BAM
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.bam.bai") \
    into GENOME_ALIGNMENT_BAI
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

    """
    hisat2 ${params.alignment_arguments.genome} \
           -x ${genome_base} -U ${fastq} \
           -p ${task.cpus} \
           --no-softclip \
           --al-gz ${sample}.${index}.genome_alignment.aligned.fastq.gz \
           --un-gz ${sample}.${index}.genome_alignment.unaligned.fastq.gz \
               2> ${sample}.${index}.genome_alignment.log \
           | samtools view -bS - \
           | samtools sort -@ ${task.cpus} -o ${sample}.${index}.genome_alignment.bam \
           && samtools index -@ ${task.cpus} ${sample}.${index}.genome_alignment.bam \
           && samtools idxstats -@ ${task.cpus} ${sample}.${index}.genome_alignment.bam > \
              ${sample}.${index}.genome_alignment.stats \
           && rfc hisat2-log-to-csv \
                  -l ${sample}.${index}.genome_alignment.log \
                  -n ${sample} -p genome \
                  -o ${sample}.${index}.genome_alignment.csv
    """
}

GENOME_ALIGNMENT_ALIGNED.into { GENOME_ALIGNMENT_ALIGNED_FASTQ_READ_LENGTH
    GENOME_ALIGNMENT_ALIGNED_MERGE
    GENOME_ALIGNMENT_ALIGNED_FASTQ_FASTQC }

GENOME_ALIGNMENT_UNALIGNED.into { GENOME_ALIGNMENT_UNALIGNED_FASTQ_READ_LENGTH
    GENOME_ALIGNMENT_UNALIGNED_MERGE
    GENOME_ALIGNMENT_UNALIGNED_FASTQ_FASTQC
    FOR_POST_GENOME }

// GENOME ALIGNMENT
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
process genome_quality_filter {
    storeDir get_storedir('quality_filter')

    input:
        set val(sample), val(index), file(bam) from GENOME_ALIGNMENT_BAM_FOR_QPASS

    output:
        set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.bam") into GENOME_ALIGNMENT_QPASS_BAM
        set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.bam.bai") into GENOME_ALIGNMENT_QPASS_BAI
        set val(sample), val(index), file("${sample}.${index}.genome.qpass.count") into GENOME_QPASS_COUNTS
        set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.stats") into GENOME_QPASS_STATS

        """
    samtools view -bq ${params.mapping_quality_cutoff} ${bam} > ${sample}.${index}.genome_alignment.qpass.bam && \
    samtools index ${sample}.${index}.genome_alignment.qpass.bam
    samtools view -c ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome.qpass.count
    samtools flagstat ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome_alignment.qpass.stats
    """
}

// Always define PSITE-ref channel to avoid undefined variable errors
if (params.containsKey('psite_offset')) {
    GENOME_ALIGNMENT_QPASS_BAM.into {
        GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE
        GENOME_ALIGNMENT_QPASS_BAM_FOR_STATS
        GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_REF
    }
} else {
    GENOME_ALIGNMENT_QPASS_BAM.into {
        GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE
        GENOME_ALIGNMENT_QPASS_BAM_FOR_STATS
    }
    Channel.empty().set { GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_REF }
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
    dedup_method == 'umi_tools' || dedup_method == 'none'

        """
  if [ "${dedup_method}" == "umi_tools" ]; then
    samtools merge ${sample}.genome.qpass.merged.bam ${bam_files} && samtools index ${sample}.genome.qpass.merged.bam
    # Create empty dummy files for the none output
    touch ${sample}.genome.qpass.merged.none.bam && touch ${sample}.genome.qpass.merged.none.bam.bai
  else
    # Create empty dummy files for the umi output
    touch ${sample}.genome.qpass.merged.bam && touch ${sample}.genome.qpass.merged.bam.bai
    samtools merge ${sample}.genome.qpass.merged.none.bam ${bam_files} && samtools index ${sample}.genome.qpass.merged.none.bam
  fi
  """
}

if (params.containsKey('psite_offset') && dedup_method == 'position') {
    GENOME_MERGED_BAM_QPASS_PRE_DEDUP_RAW.into {
        GENOME_MERGED_BAM_QPASS_PRE_DEDUP
        GENOME_MERGED_BAM_QPASS_FOR_PSITE
    }
} else {
    GENOME_MERGED_BAM_QPASS_PRE_DEDUP_RAW.set { GENOME_MERGED_BAM_QPASS_PRE_DEDUP }
}

process individual_genome_bam_to_bed {
    storeDir get_storedir('bam_to_bed') + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bam) from GENOME_ALIGNMENT_BAM_FOR_BED

    output:
        set val(sample), val(index), file("${sample}.${index}.genome.bed") into GENOME_INDIVIDUAL_BED_PRE_SPLIT
        set val(sample), val(index), file("${sample}.${index}.genome_nodedup_count.txt") \
       into GENOME_INDIVIDUAL_NODEDUP_COUNT

        """
    if [ `samtools view -c ${bam}` -eq 0 ];
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
    samtools merge ${sample}.genome.bam ${bam} && samtools index ${sample}.genome.bam && \
    zcat ${aligned_fastq} | gzip -c > ${sample}.genome.aligned.fastq.gz && \
    zcat ${unaligned_fastq} | gzip -c > ${sample}.genome.unaligned.fastq.gz && \
    rfc merge-hisat2-logs -o ${sample}.genome.log ${alignment_log} && \
    rfc hisat2-log-to-csv -l ${sample}.genome.log -n ${sample} -p genome \
                      -o ${sample}.genome.csv
    """
}

///////////////////////////////////////////////////////////////////////////////
//          G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////

// Set genome BAM channel for bed conversion and bigwig reference
GENOME_ALIGNMENT_MERGED_BAM.into {
    GENOME_BAM_FOR_DEDUP
    GENOME_ALIGNMENT_MERGED_BAM_FOR_BIGWIG_REF
}

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
    if [ `samtools view -c ${bam}` -eq 0 ];
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
    GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_UMI
    GENOME_BED_FOR_DEDUP_PRE_FOR_STATS
}

// Position-based genome deduplication
process genome_deduplicate_position {
    storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

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

if (params.containsKey('psite_offset') && dedup_method == 'position') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_PSITE
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS
    }
} else if (dedup_method == 'position') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS
    }
    // Create empty channel for PSITE when psite_offset is not configured
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_PSITE }
} else if (dedup_method == 'umi_tools') {
    // For umi_tools, the POSITION_RAW channel is not used, set empty for the mix
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX }
    // Create empty channel for PSITE when dedup_method is umi_tools
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_PSITE }
    // Create empty channel for STATS when dedup_method is umi_tools
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS }
} else {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX }
    // Create empty channel for PSITE when dedup_method is not position
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_PSITE }
    // Create empty channel for STATS when dedup_method is not position
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS }
}

process genome_deduplicate_umi_tools {
    storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam), file(bai) from GENOME_MERGED_BAM_QPASS_PRE_DEDUP

    output:
    set val(sample), file("${sample}.dedup.bam") into GENOME_UMI_TOOLS_DEDUP_BAM
    set val(sample), file("${sample}.dedup.log") into GENOME_UMI_TOOLS_DEDUP_LOG
    set val(sample), file("${sample}.dedup.stats_edit_distance.tsv"),
               file("${sample}.dedup.stats_per_umi_per_position.tsv"),
               file("${sample}.dedup.stats_per_umi.tsv") \
                    into GENOME_UMI_TOOLS_DEDUP_STATS
    set val(sample), file("${sample}.dedup.bam") into GENOME_UMI_DEDUP_BAM_FOR_BED
    set val(sample), file("${sample}.dedup.bam") into GENOME_UMI_DEDUP_BAM_FOR_SEPARATION

    when:
    dedup_method == 'umi_tools'

    """
    umi_tools dedup ${params.get('umi_tools_dedup_arguments', '')} \
        -I ${bam} --output-stats=${sample}.dedup.stats -S ${sample}.dedup.bam -L ${sample}.dedup.log
    """
}

process genome_umi_dedup_bam_to_bed {
    storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam) from GENOME_UMI_DEDUP_BAM_FOR_BED

    output:
    set val(sample), file("${sample}.dedup.bed") into GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI

    when:
    dedup_method == 'umi_tools'

    """
    bamToBed -i ${bam} > ${sample}.dedup.bed
    """
}
GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX.mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI)
.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP }

///////////////////////////////////////////////////////////////////////////////
//          P - S I T E   O F F S E T   C O R R E C T I O N
///////////////////////////////////////////////////////////////////////////////

// Initialize empty PSITE channels for workflows that don't have process outputs
Channel.empty().set { GENOME_PSITE_COUNTS_FROM_BAM_OUTPUT }
Channel.empty().set { GENOME_MERGED_PSITE_STATS_UMI }
Channel.empty().set { GENOME_MERGED_PSITE_STATS_NONE }

if (params?.psite_offset?.offset_file) {
    def offsetPath = params.psite_offset.offset_file.toString()

    // Create all channels regardless of dedup_method to avoid undefined variable errors
    Channel.fromPath(offsetPath).set { PSITE_OFFSET_FILE_UMI }
    Channel.fromPath(offsetPath).set { PSITE_OFFSET_FILE_POSITION }
    Channel.fromPath(offsetPath).set { PSITE_OFFSET_FILE_NONE }
}

Channel.empty().set { GENOME_QPASS_BAM_SINGLE_FOR_PSITE }

if (params.containsKey('psite_offset') && dedup_method == 'position') {
    GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_REF
    .map { sample, index, bam -> [sample, bam] }
    .groupTuple()
    .map { sample, bams -> [sample, bams[0]] }
    .set { GENOME_QPASS_BAM_SINGLE_FOR_PSITE }
}

// Defensive empty channels for optional features
Channel.empty().set { POST_GENOME_INDEX }
Channel.empty().set { POST_GENOME_ALIGNMENT_BAM }
Channel.empty().set { TRANSCRIPTOME_ALIGNMENT_BAM }

process apply_psite_correction_position_bed {
    storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(dedup_bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_PSITE
        file(offset_csv) from PSITE_OFFSET_FILE_POSITION.first()

    output:
        set val(sample), file("${sample}.psite.bed") into GENOME_PSITE_CORRECTED_BED_POSITION

    when:
    params.containsKey('psite_offset') && dedup_method == 'position'

    script:
    def experiment_id = params.psite_offset.experiment_mapping[sample]

    if (experiment_id == null) {
        error "No experiment_mapping found for sample '${sample}' in psite_offset configuration"
    }

    """
    python3 ${baseDir}/scripts/apply_psite_offsets_bed.py \\
        -i ${dedup_bed} \\
        -o ${sample}.psite.bed \\
        -c ${offset_csv} \\
        -e ${experiment_id} \\
        -s ${sample}
    """
}

process apply_psite_correction_umi {
    storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(dedup_bam) from GENOME_UMI_TOOLS_DEDUP_BAM
        file(offset_csv) from PSITE_OFFSET_FILE_UMI

    output:
        set val(sample), file("${sample}.psite.bam"), file("${sample}.psite.bam.bai") into GENOME_PSITE_CORRECTED_BAM_UMI
        set val(sample), file("${sample}.psite.bed") into GENOME_PSITE_CORRECTED_BED_UMI

    when:
    params.containsKey('psite_offset') && dedup_method == 'umi_tools'

    script:
    def experiment_id = params.psite_offset.experiment_mapping[sample]

    if (experiment_id == null) {
        error "No experiment_mapping found for sample '${sample}' in psite_offset configuration"
    }

    """
    python3 ${baseDir}/scripts/apply_psite_offsets.py \\
        -i ${dedup_bam} \\
        -o ${sample}.psite.bam \\
        -c ${offset_csv} \\
        -e ${experiment_id} \\
        -s ${sample} \\
        --index

    bamToBed -i ${sample}.psite.bam > ${sample}.psite.bed
    """
}

process apply_psite_correction_none {
    storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory
    input:
        set val(sample), file(qpass_bam), file(qpass_bai) from GENOME_MERGED_BAM_QPASS_NONE
        file(offset_csv) from PSITE_OFFSET_FILE_NONE.first()
    output:
        set val(sample), file("${sample}.psite.bam"), file("${sample}.psite.bam.bai") into GENOME_PSITE_CORRECTED_BAM_NONE
        set val(sample), file("${sample}.psite.bed") into GENOME_PSITE_CORRECTED_BED_NONE
    when:
    params.containsKey('psite_offset') && dedup_method == 'none'
    script:
    def experiment_id = params.psite_offset.experiment_mapping[sample]
    if (experiment_id == null) {
        error "No experiment_mapping found for sample '${sample}' in psite_offset configuration"
    }
    """
    python3 ${baseDir}/scripts/apply_psite_offsets.py \\
        -i ${qpass_bam} \\
        -o ${sample}.psite.bam \\
        -c ${offset_csv} \\
        -e ${experiment_id} \\
        -s ${sample} \\
        --index

    bamToBed -i ${sample}.psite.bam > ${sample}.psite.bed
    """
}


if (params.containsKey('psite_offset') && dedup_method == 'umi_tools') {
    GENOME_PSITE_CORRECTED_BAM_UMI.into {
        GENOME_BAM_FOR_DOWNSTREAM
        GENOME_BAM_FOR_BIGWIG_FINAL
        GENOME_PSITE_COUNTS_FROM_BAM
    }
    // Set up PSITE channels for UMI-tools workflow with P-site correction
    GENOME_PSITE_CORRECTED_BED_UMI.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL
        GENOME_PSITE_BED_FOR_BIGWIG
        GENOME_PSITE_BED_FOR_STATS
    }

    // Generate P-site counts from psite.bam files for UMI-tools deduplication
    process generate_psite_counts_from_bam {
        storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(psite_bam), file(psite_bai) from GENOME_PSITE_COUNTS_FROM_BAM

        output:
            set val(sample), file("${sample}.psite_count.txt") into GENOME_PSITE_COUNTS_FROM_BAM_OUTPUT

        """
        samtools view -c ${psite_bam} > ${sample}.psite_count.txt
        """
    }

    
    } else if (params.containsKey('psite_offset') && dedup_method == 'position') {
    GENOME_PSITE_CORRECTED_BED_POSITION.into {
        GENOME_PSITE_BED_FOR_BIGWIG
        GENOME_PSITE_BED_FOR_STATS
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL
    }

  
    // Generate post-dedup counts from BED files for position deduplication (UMI-tools approach)
    process generate_post_dedup_counts_from_bed_position {
        storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(post_dedup_bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS

        output:
        set val(sample), file("${sample}.count_after_dedup.txt") into GENOME_POSITION_DEDUP_COUNTS

        """
        wc -l < ${post_dedup_bed} > ${sample}.count_after_dedup.txt
        """
    }

    // Split merged post-dedup counts into individual format for position deduplication
    process split_post_dedup_counts_for_position {

        input:
        set val(sample), file(count_file) from GENOME_POSITION_DEDUP_COUNTS

        output:
        set val(sample), val(1), file("${sample}.1.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING

        """
        # For position dedup, we use the same count for all indices since dedup is done after merging
        # Get the count from the merged file
        count=\$(cat ${count_file})

        # Create individual count file with index 1
        echo "\${count}" > ${sample}.1.count_after_dedup.txt
        """
    }

    // Generate P-site counts from BED files for position deduplication
    process generate_psite_counts_from_bed_position {
        storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(psite_bed) from GENOME_PSITE_BED_FOR_STATS

        output:
        set val(sample), file("${sample}.psite_count.txt") into GENOME_PSITE_COUNTS_FROM_BED_POSITION_OUTPUT

        """
        wc -l < ${psite_bed} > ${sample}.psite_count.txt
        """
    }

Channel.empty().set { GENOME_BAM_FOR_BIGWIG_FINAL }
} else if (params.containsKey('psite_offset') && dedup_method == 'none') {
    GENOME_PSITE_CORRECTED_BAM_NONE.into {
        GENOME_BAM_FOR_BIGWIG_FINAL
        GENOME_PSITE_COUNTS_FROM_BAM
    }
    GENOME_PSITE_CORRECTED_BED_NONE.into {
        GENOME_PSITE_BED_FOR_BIGWIG
        GENOME_PSITE_BED_FOR_STATS
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL
    }

    // Generate P-site counts from psite.bam files for none deduplication
    process generate_psite_counts_from_bam_none {
        storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(psite_bam), file(psite_bai) from GENOME_PSITE_COUNTS_FROM_BAM

        output:
            set val(sample), file("${sample}.psite_count.txt") into GENOME_PSITE_COUNTS_FROM_BAM_OUTPUT_NONE

        """
        samtools view -c ${psite_bam} > ${sample}.psite_count.txt
        """
    }

    // Generate P-site statistics for merged samples (none dedup with P-site)
    // Old process removed - P-site stats are now handled by add_psite_stats_to_merged process
} else if (dedup_method == 'umi_tools') {
    GENOME_UMI_TOOLS_DEDUP_BAM.into {
        GENOME_BAM_FOR_DOWNSTREAM
        GENOME_DEDUP_BAM_FOR_INDEXING
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL }

    process index_dedup_bam_for_bigwig {
        storeDir get_storedir('deduplication') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(bam) from GENOME_DEDUP_BAM_FOR_INDEXING

        output:
            set val(sample), file(bam), file("${bam}.bai") into GENOME_BAM_FOR_BIGWIG_FINAL

            """
        samtools index -@ ${task.cpus} ${bam}
        """
    }
} else {
    // Position-based without P-site correction
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL }
    // Create empty BAM channel since position-based uses BED workflow
    Channel.empty().set { GENOME_BAM_FOR_BIGWIG_FINAL }
}

// Only split BAM to individual files for UMI-tools workflow
if (dedup_method == 'umi_tools') {
    GENOME_BAM_FOR_DOWNSTREAM
    .cross(GENOME_INDIVIDUAL_BED_FOR_SPLITTING.map { sample, index, bed -> [sample, index] }.unique())
    .map { merged_data, individual_data ->
        [individual_data[0], individual_data[1], merged_data[1]]
    }
    .set { GENOME_DEDUP_BAM_FOR_SPLITTING }

    process split_genome_dedup_bam_to_individual {
        storeDir get_storedir('deduplication') + '/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(merged_bam) from GENOME_DEDUP_BAM_FOR_SPLITTING

        output:
            set val(sample), val(index), file("${sample}.${index}.post_dedup.bam"), \
                             file("${sample}.${index}.post_dedup.bam.bai") into GENOME_INDIVIDUAL_DEDUP_BAM
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING

        """
        samtools view -B -r ${sample}.${index} ${merged_bam} -o ${sample}.${index}.post_dedup.bam
        samtools view -c ${sample}.${index}.post_dedup.bam > ${sample}.${index}.count_after_dedup.txt
        samtools index -@ ${task.cpus} ${sample}.${index}.post_dedup.bam
        """
    }
}

if (dedup_method == 'umi_tools') {
    // Use individual dedup counts from BAM splitting (has sample, index, count structure)
    GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING.set { GENOME_INDIVIDUAL_DEDUP_COUNT }

    GENOME_INDIVIDUAL_DEDUP_BAM.set { GENOME_INDIVIDUAL_DEDUP_BAM_FOR_STATS }
}
else if (dedup_method == 'position') {
    // For position-based dedup, post-dedup counts are generated by generate_post_dedup_counts_from_bed_position process
    // This follows the same approach as UMI-tools deduplication
    GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
}
else {
    // No deduplication - use individual nodedup counts from BAM to BED
    GENOME_INDIVIDUAL_NODEDUP_COUNT
    .set { GENOME_INDIVIDUAL_DEDUP_COUNT }
}

///////////////////////////////////////////////////////////////////////////////
//  R I B O - S E Q   B I G W I G   G E N E R A T I O N
///////////////////////////////////////////////////////////////////////////////

// BAM-based bigWig generation (for UMI-tools with P-site correction)
process genome_create_strand_specific_bigwigs {
    storeDir get_storedir('deduplication') + '/bigwigs/' + params.output.merged_lane_directory

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_bigwig'

    input:
    set val(sample), file(bam), file(bai) from GENOME_BAM_FOR_BIGWIG_FINAL

    output:
    set val(sample), file("${sample}.psite.plus.bigWig"), \
                     file("${sample}.psite.minus.bigWig") \
        into GENOME_STRAND_SPECIFIC_BIGWIGS

    when:
    (dedup_method == 'umi_tools' || dedup_method == 'none') && params.containsKey('psite_offset')

    """
    bamCoverage -b ${bam} -o ${sample}.psite.plus.bigWig \
        --filterRNAstrand reverse --binSize 1 -p ${task.cpus}

    bamCoverage -b ${bam} -o ${sample}.psite.minus.bigWig \
        --filterRNAstrand forward --binSize 1 -p ${task.cpus}
    """
}

// BED-based bigWig generation (for position-based dedup with P-site correction)
process genome_create_psite_bigwigs_from_bed {
    storeDir get_storedir('deduplication') + '/bigwigs/' + params.output.merged_lane_directory

    input:
    set val(sample), file(psite_bed) from GENOME_PSITE_BED_FOR_BIGWIG
    set val(sample_bam), file(ref_bam) from GENOME_ALIGNMENT_MERGED_BAM_FOR_BIGWIG_REF.first()

    output:
    set val(sample), file("${sample}.psite.plus.bigWig"), \
                     file("${sample}.psite.minus.bigWig") \
        into GENOME_PSITE_BIGWIGS_FROM_BED

    when:
    dedup_method == 'position' && params.containsKey('psite_offset')

    """
    # Extract genome sizes from BAM header
    samtools view -H ${ref_bam} | grep '@SQ' | awk '{print substr(\$2,4)"\\t"substr(\$3,4)}' > genome.sizes

    # Separate by strand and create bedGraph files
    awk '\$6=="+"' ${psite_bed} | sort -k1,1 -k2,2n > plus.sorted.bed
    awk '\$6=="-"' ${psite_bed} | sort -k1,1 -k2,2n > minus.sorted.bed

    # Create coverage bedGraph (using genomecov with -bg for bedGraph output)
    bedtools genomecov -i plus.sorted.bed -g genome.sizes -bg > plus.bedGraph
    bedtools genomecov -i minus.sorted.bed -g genome.sizes -bg > minus.bedGraph

    # Convert bedGraph to bigWig
    bedGraphToBigWig plus.bedGraph genome.sizes ${sample}.psite.plus.bigWig
    bedGraphToBigWig minus.bedGraph genome.sizes ${sample}.psite.minus.bigWig
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

// Create indexed channels for clip and filter logs for genome stats
CLIP_LOG.map { sample, index, clip_log -> [ [sample, index], clip_log ] }
        .set { CLIP_LOG_INDEXED_FOR_GENOME }
FILTER_LOG.map { sample, index, filter_log -> [ [sample, index], filter_log ] }
          .set { FILTER_LOG_INDEXED_FOR_GENOME }

///////////////////////////////////////////////////////////////////////////////
/* GENOME STATS */

// Create indexed channels for genome stats
GENOME_ALIGNMENT_LOG_STATS
  .map { sample, index, log_file -> [ [sample, index], log_file] }
  .set { GENOME_ALIGNMENT_LOG_STATS_INDEXED }

GENOME_QPASS_COUNTS
  .map { sample, index, count_file -> [ [sample, index], count_file] }
  .set { GENOME_QPASS_COUNTS_INDEXED }

GENOME_INDIVIDUAL_DEDUP_COUNT
  .map { sample, index, count_file -> [ [sample, index], count_file] }
  .set { GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED }

// Create genome alignment stats input channel
CLIP_LOG_INDEXED_FOR_GENOME.join(FILTER_LOG_INDEXED_FOR_GENOME)
            .join(GENOME_ALIGNMENT_LOG_STATS_INDEXED)
            .join(GENOME_QPASS_COUNTS_INDEXED)
            .join(GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED)
            .flatten()
            .collate(7)
            .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }

// INDIVIDUAL GENOME ALIGNMENT STATS
process individual_genome_alignment_stats {
    storeDir get_storedir('stats') + '/genome/individual'

input:
    set val(sample), val(index), file(clip_log), file(filter_log),\
    file(genome_log), file(qpass_count),\
    file(dedup_count)\
    from GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT

output:
    set val(sample), val(index), file("${sample}.${index}.genome_individual.csv") \
   into GENOME_INDIVIDUAL_ALIGNMENT_STATS

    """
    # Build the command
    cmd="rfc compile-step-stats \\
      -n ${sample}.${index} \\
      -c ${clip_log} \\
      -f ${filter_log} \\
      -g ${genome_log} \\
      -q ${qpass_count} \\
      -d ${dedup_count} \\
      -l genome \\
      -o ${sample}.${index}.genome_individual.csv"
    if [[ -n "${params.psite_offset}" ]]; then
        psite_count_path="${get_storedir('deduplication')}/${params.output.merged_lane_directory}/${sample}.psite_count.txt"
        if [ -f "\$psite_count_path" ]; then
            cmd="\${cmd} -p \$psite_count_path"
            echo "Using P-site count: \$psite_count_path"
        else
            echo "P-site count file not found: \$psite_count_path"
        fi
    fi

    eval "\$cmd"
    """
}

// Split genome individual stats
GENOME_INDIVIDUAL_ALIGNMENT_STATS
  .into { GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
      GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING}

// Collect genome individual stats
GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
  .map { sample, index, stats_file -> stats_file }
  .toSortedList().set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED }

// COLLECT INDIVIDUAL GENOME ALIGNMENT STATS
process combine_individual_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/individual'

input:
    file(stat_table) from GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

output:
    file('genome_individual_essential.csv') \
      into COMBINED_INDIVIDUAL_GENOME_ALIGNMENT_STATS

    script:
    if (stat_table.size() == 0) {
        // Create empty stats file when no input is available
        """
        echo "No individual statistics data available" > genome_individual_essential.csv
        echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> genome_individual_essential.csv
        """
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

// MERGED GENOME ALIGNMENT STATS
GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING
  .map { sample, index, file -> [ sample, file ] }
  .groupTuple()
  .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED }

process sum_individual_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/merged'

input:
    set val(sample), file(stat_files) from GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED

output:
    set val(sample), file("${sample}.genome_merged.csv") into GENOME_MERGED_ALIGNMENT_STATS

    """
    rfc sum-stats -n ${sample}\
  -o ${sample}.genome_merged.csv ${stat_files}
    """
}

GENOME_MERGED_ALIGNMENT_STATS.set { GENOME_MERGED_ALIGNMENT_STATS_FINAL }

// Split and collect merged genome stats
GENOME_MERGED_ALIGNMENT_STATS_FINAL
  .into { GENOME_MERGED_ALIGNMENT_STATS_FOR_COLLECTION
      GENOME_MERGED_ALIGNMENT_STATS_FOR_GROUPING}

GENOME_MERGED_ALIGNMENT_STATS_FOR_COLLECTION
  .map { sample, stats_file -> stats_file }
  .toSortedList().set { GENOME_MERGED_ALIGNMENT_STATS_COLLECTED }

// COMBINE MERGED GENOME ALIGNMENT STATS
process combine_merged_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/merged'

input:
    file(stat_files) from GENOME_MERGED_ALIGNMENT_STATS_COLLECTED

output:
    file('genome_merged_essential.csv') into COMBINED_MERGED_GENOME_ALIGNMENT_STATS

    script:
    if (stat_files.size() == 0) {
        // Create empty stats file when no input is available
        """
        echo "No merged genome statistics data available" > genome_merged_essential.csv
        echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup,genome_after_psite" >> genome_merged_essential.csv
        """
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

///////////////////////////////////////////////////////////////////////////////
/* METADATA CHANNELS */
do_metadata = params.get('do_metadata', false) && params.get('input', [:]).get('metadata', false)

if (do_metadata) {
    meta_base = params.input.metadata.get('base', '')
    if (meta_base != '' && !meta_base.endsWith('/')) {
        meta_base = "${meta_base}/"
    }

    Channel.from(params.input.metadata.files.collect { k, v ->
                           [k, file("${meta_base}${v}") ] })
                                             .into { METADATA_PRE; METADATA_PRE_VERBOSE }

    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL
                   .join(METADATA_PRE, remainder: true)
                   .into { METADATA_RIBO; METADATA_VERBOSE }
}
else {
    METADATA_PRE = Channel.from([null])
    METADATA_PRE_VERBOSE = Channel.from([null])

    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL
  .map { sample, bed -> [sample, bed, null] }
  .into { METADATA_RIBO; METADATA_VERBOSE }
}

if (params.get('input', [:]).get('root_meta', false)) {
    ROOT_META = Channel.from([file(params.input.root_meta)])
}
else {
    ROOT_META = Channel.from([null])
}

// METADATA CHANNELS
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/* CREATE RIBO FILES */

do_create_ribo = params.get('do_ribo_creation', true)

if (do_create_ribo) {
    Channel.from(file(params.input.reference.transcript_lengths)).
into { T_LENGTHS_FOR_RIBO; T_LENGTHS_FOR_IND_RIBO }

    Channel.from(file(params.input.reference.regions)).
into { ANNOTATION_FOR_RIBO; ANNOTATION_FOR_IND_RIBO }

    if (params.ribo.coverage) {
        coverage_argument = ''
    }
else {
        coverage_argument = '--nocoverage'
}

    process create_ribo {
        publishDir get_publishdir('ribo') + '/experiments', mode:'copy'
        storeDir get_storedir('ribo') + '/experiments'

    input:
        set val(sample), file(bed_file), file(meta_file) from METADATA_RIBO
        file(transcript_length_file) from T_LENGTHS_FOR_RIBO.first()
        file(annotation_file) from ANNOTATION_FOR_RIBO.first()
        file(root_meta_file) from ROOT_META.first()

    output:
        set val(sample), file("${sample}.ribo") into RIBO_MAIN

  script:
        if (meta_file != null) {
            sample_meta_argument = "--expmeta ${meta_file}"
        }
  else {
            sample_meta_argument = ''
  }

        if (root_meta_file == null) {
            root_meta_argument = ''
        }
  else {
            root_meta_argument = "--ribometa ${root_meta_file}"
  }

        """
    ribopy create -n ${sample} \
                 --reference ${params.ribo.ref_name} \
                 --lengths ${transcript_length_file} \
                             --annotation ${annotation_file} \
                             --radius ${params.ribo.metagene_radius} \
                             -l ${params.ribo.left_span} -r ${params.ribo.right_span} \
                             --lengthmin ${params.ribo.read_length.min} \
                             --lengthmax ${params.ribo.read_length.max} \
               ${sample_meta_argument} \
               ${root_meta_argument} \
                             ${coverage_argument} \
                             -n ${task.cpus} \
               --alignmentfile ${bed_file} \
                ${sample}.ribo
    """
    }

    RIBO_MAIN.into { RIBO_FOR_RNASEQ; RIBO_AFTER_CREATION }
} else {
    RIBO_FOR_RNASEQ = Channel.empty()
    RIBO_AFTER_CREATION = Channel.empty()
} // end if(do_create_ribo)

// CREATE RIBO FILES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/* Post Genome */

do_post_genome = params.input.reference.get('post_genome', false)

if (do_post_genome) {
    POST_GENOME_INDEX = Channel.from([[
                  params.input.reference.post_genome
                  .split('/')[-1]
                  .replaceAll('\\*$', '')
                  .replaceAll('\\.$', ''),
               file(params.input.reference.post_genome),
              ]])

    process post_genome_alignment {
        storeDir get_storedir('post_genome_alignment') + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(fastq) from FOR_POST_GENOME
        set val(post_genome_base), file(post_genome_files) from POST_GENOME_INDEX.first()

  output:
        set val(sample), val(index), file("${sample}.${index}.postgenome_alignment.bam") \
      into POST_GENOME_ALIGNMENT_BAM
        set val(sample), val(index), file("${sample}.${index}.postgenome_alignment.bam.bai") \
      into POST_GENOMEE_ALIGNMENT_BAI
        set val(sample), val(index), file("${sample}.${index}.aligned.postgenome_alignment.fastq.gz") \
      into POST_GENOME_ALIGNMENT_ALIGNED
        set val(sample), val(index), file("${sample}.${index}.unaligned.postgenome_alignment.fastq.gz") \
      into POST_GENOME_ALIGNMENT_UNALIGNED
        set val(sample), val(index), file("${sample}.${index}.postgenome_alignment.log") \
      into POST_GENOME_ALIGNMENT_LOG
        set val(sample), val(index), file("${sample}.${index}.postgenome_alignment.csv") \
      into POST_GENOME_ALIGNMENT_CSV
        set val(sample), val(index), file("${sample}.${index}.postgenome_alignment.stats") \
      into POST_GENOME_ALIGNMENT_STATS

        """
  bowtie2 ${params.alignment_arguments.transcriptome} \
          -x ${post_genome_base} -q ${fastq} \
          --threads ${task.cpus} \
          --al-gz ${sample}.${index}.aligned.postgenome_alignment.fastq.gz \
          --un-gz ${sample}.${index}.unaligned.postgenome_alignment.fastq.gz \
                     2> ${sample}.${index}.postgenome_alignment.log \
          | samtools view -bS - \
          | samtools sort -@ ${task.cpus} -o ${sample}.${index}.postgenome_alignment.bam \
          && samtools index -@ ${task.cpus} ${sample}.${index}.postgenome_alignment.bam \
          && samtools idxstats -@ ${task.cpus} ${sample}.${index}.postgenome_alignment.bam  > \
             ${sample}.${index}.postgenome_alignment.stats \
          && rfc bt2-log-to-csv -o ${sample}.${index}.postgenome_alignment.csv \
                -n ${sample} -p post_genome -l ${sample}.${index}.postgenome_alignment.log
  """
    }

    POST_GENOME_ALIGNMENT_ALIGNED.into { POST_GENOME_ALIGNMENT_ALIGNED_FASTQ_READ_LENGTH
        POST_GENOME_ALIGNMENT_ALIGNED_MERGE
        POST_GENOME_ALIGNMENT_ALIGNED_FASTQ_FASTQC }

    POST_GENOME_ALIGNMENT_UNALIGNED.into { POST_GENOME_ALIGNMENT_UNALIGNED_FASTQ_READ_LENGTH
        POST_GENOME_ALIGNMENT_UNALIGNED_MERGE
        POST_GENOME_ALIGNMENT_UNALIGNED_FASTQ_FASTQC}

// POST_GENOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* MERGE POST_GENOME ALIGNMENT */
    POST_GENOME_ALIGNMENT_LOG.into { POST_GENOME_ALIGNMENT_LOG_MERGE; POST_GENOME_ALIGNMENT_LOG_TABLE }

    POST_GENOME_ALIGNMENT_BAM.map { sample, index, bam -> [sample, bam] }.groupTuple()
    .set { POST_GENOME_ALIGNMENT_GROUPED_BAM }

    POST_GENOME_ALIGNMENT_ALIGNED_MERGE.map { sample, index, fastq -> [sample, fastq] }.groupTuple()
    .set { POST_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ }

    POST_GENOME_ALIGNMENT_UNALIGNED_MERGE.map { sample, index, fastq -> [sample, fastq] }.groupTuple()
    .set { POST_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ }
    POST_GENOME_ALIGNMENT_LOG_MERGE.map { sample, index, log -> [sample, log] }.groupTuple()
    .set { POST_GENOME_ALIGNMENT_GROUPED_LOG }

    POST_GENOME_ALIGNMENT_GROUPED_BAM.join(POST_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ)
                            .join(POST_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ)
                            .join(POST_GENOME_ALIGNMENT_GROUPED_LOG)
                            .set { POST_GENOME_ALIGNMENT_GROUPED_JOINT }

    process merge_post_genome_alignment {
        storeDir get_storedir('post_genome_alignment') + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(bam), file(aligned_fastq), \
          file(unaligned_fastq), file(alignment_log) from POST_GENOME_ALIGNMENT_GROUPED_JOINT

    output:
        set val(sample), file("${sample}.post_genome.bam") \
      into POST_GENOME_ALIGNMENT_MERGED_BAM
        set val(sample), file("${sample}.post_genome.bam.bai") \
      into POST_GENOME_ALIGNMENT_MERGED_BAI
        set val(sample), file("${sample}.post_genome.aligned.fastq.gz") \
      into POST_GENOME_ALIGNMENT_MERGED_ALIGNED_FASTQ
        set val(sample), file("${sample}.post_genome.unaligned.fastq.gz") \
      into POST_GENOME_ALIGNMENT_MERGED_UNALIGNED_FASTQ
        set val(sample), file("${sample}.post_genome.log") \
      into POST_GENOME_ALIGNMENT_MERGED_LOG
        set val(sample), file("${sample}.post_genome.csv") \
      into POST_GENOME_ALIGNMENT_MERGED_CSV

        """
    samtools merge ${sample}.post_genome.bam ${bam} && samtools index ${sample}.post_genome.bam && \
    zcat ${aligned_fastq} | gzip -c > ${sample}.post_genome.aligned.fastq.gz && \
    zcat ${unaligned_fastq} | gzip -c > ${sample}.post_genome.unaligned.fastq.gz && \
    rfc merge bowtie2-logs -o ${sample}.post_genome.log ${alignment_log} && \
    rfc bt2-log-to-csv -o ${sample}.post_genome.csv \
          -n ${sample} -p post_genome -l ${sample}.post_genome.log
    """
    }

    POST_GENOME_ALIGNMENT_CSV
   .map { sample, index, stats_file -> stats_file }
   .toSortedList().set { POST_GENOME_ALIGNMENT_CSV_INDIVIDUAL_LIST }

    POST_GENOME_ALIGNMENT_MERGED_CSV
   .map { sample, stats_file -> stats_file }
   .toSortedList().set { POST_GENOME_ALIGNMENT_CSV_MERGED_LIST }

    process combine_individual_postgenome_stats {
        storeDir get_storedir('post_genome_alignment') + '/logs'

  input:
        file(stats_input_files) from POST_GENOME_ALIGNMENT_CSV_INDIVIDUAL_LIST
        file(stats_input_files_merged) from POST_GENOME_ALIGNMENT_CSV_MERGED_LIST

  output:
        file('postgenome_individual_stats.csv') \
        into POST_GENOME_ALIGNMENT_CSV_INDIVIDUAL_COMBINED
        file('postgenome_merged_stats.csv') \
        into POST_GENOME_ALIGNMENT_CSV_MERGED_COMBINED

        """
  rfc merge overall-stats -o postgenome_individual_stats.csv ${stats_input_files} ; \
  rfc merge overall-stats -o postgenome_merged_stats.csv ${stats_input_files_merged}
  """
    }

// MERGE POST GENOME ALIGNMENT
////////////////////////////////////////////////////////////////////////////////
} 

// Genome is now the only alignment method - use genome stats as final
COMBINED_INDIVIDUAL_GENOME_ALIGNMENT_STATS.set { FINAL_INDIVIDUAL_STATS }

// Final stats channels are set above in the genome/transcriptome conditional blocks

// Use final stats directly for publishing
FINAL_INDIVIDUAL_STATS.set { FINAL_INDIVIDUAL_STATS_FOR_PUBLISH }


process publish_stats {
    publishDir get_publishdir('stats'), mode: 'copy'

    executor 'local'

  input:
    file(individual_stats) from FINAL_INDIVIDUAL_STATS_FOR_PUBLISH

  output:
    file('individual_stats.csv') into INDIVIDUAL_STATS_PUBLISHED

    """
  cp ${individual_stats} individual_stats.csv
  # Merged sample files are no longer generated - using individual files only
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

    Channel.from(params.get('rnaseq', [:]).fastq.collect { k, v ->
        v.collect { z -> [k, v.indexOf(z) + 1,
                                                   file("${rnaseq_fastq_base}${z}")] }  })
    .flatten().collate(3).into {  RNASEQ_FASTQ
                             RNASEQ_FASTQ_VERBOSE
                             RNASEQ_FASTQ_FASTQC
                             RNASEQ_FASTQ_CLIP
                             RNASEQ_FASTQ_EXISTENCE}

    if (params.do_check_file_existence) {
        // Make Sure Fastq Files Exist
        RNASEQ_FASTQ_EXISTENCE
  .map { sample, index, this_file -> file_exists(this_file) }
    }

    process rnaseq_raw_fastqc {
        publishDir get_publishdir('fastqc', true), mode: 'copy'

      input:
        set val(sample), val(index), file(fastq) from RNASEQ_FASTQ_FASTQC

      output:
        set val(sample), file("${sample}.${index}_fastqc.html"),
        file("${sample}.${index}_fastqc.zip") into RNASEQ_FASTQC_OUT

    when:
    params.do_fastqc && do_rnaseq

        """
      if [ ! -f ${sample}.${index}.fastq.gz ]; then
         ln -s $fastq ${sample}.${index}.fastq.gz
      fi
      fastqc ${sample}.${index}.fastq.gz --outdir=\$PWD -t ${task.cpus}
      """
    }

    process rnaseq_clip {
        storeDir get_storedir('clip', true)

  input:
        set val(sample), val(index), file(fastq) from RNASEQ_FASTQ_CLIP

  output:
        set val(sample), val(index), file("${sample}.${index}.clipped.fastq.gz") \
                                                      into RNASEQ_CLIP_OUT
        set val(sample), val(index), file("${sample}.${index}.clipped.log") \
                                                      into RNASEQ_CLIP_LOG

        """
  cutadapt --cores=${task.cpus} ${params.rnaseq.clip_arguments} ${fastq} 2>${sample}.${index}.clipped.log  \
   | gzip -c  > ${sample}.${index}.clipped.fastq.gz
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
        set val(sample), val(index), file(fastq) \
                          from RNASEQ_CLIP_OUT
        set val(bowtie2_index_base), file(bowtie2_index_files) \
                          from RNASEQ_FILTER_INDEX.first()

    output:
        set val(sample), val(index), file("${sample}.${index}.filter.bam") \
        into RNASEQ_FILTER_BAM
        set val(sample), val(index), file("${sample}.${index}.filter.bam.bai") \
        into RNASEQ_FILTER_BAI
        set val(sample), val(index), file("${sample}.${index}.aligned.filter.fastq.gz") \
        into RNASEQ_FILTER_ALIGNED
        set val(sample), val(index), file("${sample}.${index}.unaligned.filter.fastq.gz") \
        into RNASEQ_FILTER_UNALIGNED
        set val(sample), val(index), file("${sample}.${index}.filter.log") \
        into RNASEQ_FILTER_LOG
        set val(sample), val(index), file("${sample}.${index}.filter.stats") \
        into RNASEQ_FILTER_STATS

        """
    bowtie2 ${params.rnaseq.filter_arguments} \
            -x ${bowtie2_index_base} -q ${fastq} \
            --threads ${task.cpus} \
            --al-gz ${sample}.${index}.aligned.filter.fastq.gz \
            --un-gz ${sample}.${index}.unaligned.filter.fastq.gz \
                     2> ${sample}.${index}.filter.log \
            | samtools view -bS - \
            | samtools sort -@ ${task.cpus} -o ${sample}.${index}.filter.bam \
            && samtools index -@ ${task.cpus} ${sample}.${index}.filter.bam \
            && samtools idxstats -@ ${task.cpus} ${sample}.${index}.filter.bam  > \
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

// RNA-seq genome alignment (genome is now always enabled)
    // Check if RNA-seq genome read trimming is enabled (default: enabled)
    do_rnaseq_genome_trim = params.rnaseq.containsKey('genome_read_trim') &&
                           params.rnaseq.genome_read_trim.get('enabled', true) &&
                           params.rnaseq.genome_read_trim.containsKey('length')

    if (do_rnaseq_genome_trim) {
        // Add trimming step before genome alignment
        process rnaseq_genome_read_trim {
            // Keep trimmed files only in work directory - not needed as permanent output

        input:
            set val(sample), val(index), file(fastq) from RNASEQ_FILTER_UNALIGNED_GENOME

        output:
            set val(sample), val(index), file("${sample}.${index}.genome_trimmed.fastq.gz") into RNASEQ_GENOME_TRIMMED

        script:
            trim_length = params.rnaseq.genome_read_trim.length
            """
        cutadapt --length ${trim_length} -o ${sample}.${index}.genome_trimmed.fastq.gz ${fastq}
        """
        }

        // Use trimmed reads for genome alignment
        RNASEQ_GENOME_INPUT_CHANNEL = RNASEQ_GENOME_TRIMMED
} else {
        // Use original rRNA-filtered reads without trimming
        RNASEQ_GENOME_INPUT_CHANNEL = RNASEQ_FILTER_UNALIGNED_GENOME
    }

    // Create separate genome index channel for RNA-seq to avoid conflict with ribo-seq
    RNASEQ_GENOME_INDEX = Channel.from([[
            params.input.reference.genome
               .split('/')[-1]
               .replaceAll('\\*$', '')
               .replaceAll('\\.$', ''),
            file(params.input.reference.genome),
           ]])

    process rnaseq_genome_alignment {
        storeDir get_storedir('genome_alignment', true) + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(fastq) from RNASEQ_GENOME_INPUT_CHANNEL
        set val(genome_base), file(genome_files) from RNASEQ_GENOME_INDEX.first()

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

        """
    hisat2 ${params.get('rnaseq', [:]).get('hisat2_arguments', params.alignment_arguments.genome)} \\
           -x ${genome_base} -U ${fastq} \\
           -p ${task.cpus} \\
           --no-softclip \\
           --al-gz ${sample}.${index}.rnaseq_genome_alignment.aligned.fastq.gz \\
           --un-gz ${sample}.${index}.rnaseq_genome_alignment.unaligned.fastq.gz \\
                2> ${sample}.${index}.rnaseq_genome_alignment.log \\
           | samtools view -bS - \\
           | samtools sort -@ ${task.cpus} -o ${sample}.${index}.rnaseq_genome_alignment.bam \\
           && samtools index -@ ${task.cpus} ${sample}.${index}.rnaseq_genome_alignment.bam \\
           && samtools idxstats -@ ${task.cpus} ${sample}.${index}.rnaseq_genome_alignment.bam > \\
              ${sample}.${index}.rnaseq_genome_alignment.stats \\
           && rfc hisat2-log-to-csv \\
                  -l ${sample}.${index}.rnaseq_genome_alignment.log \\
                  -n ${sample} -p rnaseq_genome \\
                  -o ${sample}.${index}.rnaseq_genome_alignment.csv
    """
    }

    // Split RNA-seq genome channels for different uses
    RNASEQ_GENOME_ALIGNMENT_BAM.into {
        RNASEQ_GENOME_ALIGNMENT_BAM_FOR_QUALITY
        RNASEQ_GENOME_ALIGNMENT_BAM_FOR_MERGE
    }

    RNASEQ_GENOME_ALIGNMENT_LOG.into {
        RNASEQ_GENOME_ALIGNMENT_LOG_FOR_MERGE
        RNASEQ_GENOME_ALIGNMENT_LOG_FOR_STATS
    }

    // RNA-seq genome quality filtering (following ribo-seq genome pattern)
    process rnaseq_genome_quality_filter {
        storeDir get_storedir('quality_filter', true)

    input:
        set val(sample), val(index), file(bam) from RNASEQ_GENOME_ALIGNMENT_BAM_FOR_QUALITY

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.qpass.bam") into RNASEQ_GENOME_ALIGNMENT_QPASS_BAM
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.qpass.bam.bai") into RNASEQ_GENOME_ALIGNMENT_QPASS_BAI
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.qpass.count") into RNASEQ_GENOME_QPASS_COUNTS
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_alignment.qpass.stats") into RNASEQ_GENOME_QPASS_STATS

        """
    samtools view -bq ${params.mapping_quality_cutoff} ${bam} > ${sample}.${index}.rnaseq_genome_alignment.qpass.bam && \\
    samtools index ${sample}.${index}.rnaseq_genome_alignment.qpass.bam
    samtools view -c ${sample}.${index}.rnaseq_genome_alignment.qpass.bam > ${sample}.${index}.rnaseq_genome.qpass.count
    samtools flagstat ${sample}.${index}.rnaseq_genome_alignment.qpass.bam > ${sample}.${index}.rnaseq_genome_alignment.qpass.stats
    """
    }

    // Convert individual RNA-seq genome BAMs to BED (following ribo-seq genome pattern)
    process rnaseq_individual_genome_bam_to_bed {
        storeDir get_storedir('bam_to_bed', true) + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bam) from RNASEQ_GENOME_ALIGNMENT_QPASS_BAM

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.bed") into RNASEQ_INDIVIDUAL_GENOME_BED
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_nodedup_count.txt") into RNASEQ_INDIVIDUAL_GENOME_NODEDUP_COUNT

        """
    if [ `samtools view -c ${bam}` -eq 0 ];
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

    // Add sample index column to RNA-seq genome bed files (matching transcriptome pattern)
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

    // Group genome BED files with sample index for merging (following transcriptome pattern)
    RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED.groupTuple()
    .set { RNASEQ_GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED }

    // Merge RNA-seq genome BED files (following transcriptome pattern)
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

    // Create genome nodedup count files from individual BED files (matching transcriptome pattern)
    process rnaseq_genome_create_nodedup_counts {
    input:
        set val(sample), val(index), file(bed) from RNASEQ_GENOME_BED_FOR_STATS

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_nodedup_count.txt") into RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP

        """
    wc -l ${bed} > ${sample}.${index}.rnaseq_genome_nodedup_count.txt
    """
    }

    // RNA-seq genome individual statistics compilation (following ribo-seq genome pattern)
    // Create separate channels for RNA-seq genome statistics to avoid consumption conflicts
    RNASEQ_CLIP_LOG_INDEXED_FOR_GENOME = RNASEQ_CLIP_LOG_FOR_GENOME.map { sample, index, clip_log -> [ [sample, index], clip_log ] }
    RNASEQ_FILTER_LOG_INDEXED_FOR_GENOME = RNASEQ_FILTER_LOG_FOR_GENOME.map { sample, index, filter_log -> [ [sample, index], filter_log ] }
    RNASEQ_GENOME_ALIGNMENT_LOG_FOR_STATS.map { sample, index, genome_log -> [ [sample, index], genome_log ] }
    .set { RNASEQ_GENOME_ALIGNMENT_LOG_INDEXED }
    RNASEQ_GENOME_QPASS_COUNTS.map { sample, index, qpass_count -> [ [sample, index], qpass_count ] }
    .set { RNASEQ_GENOME_QPASS_COUNTS_INDEXED }
    RNASEQ_INDIVIDUAL_GENOME_NODEDUP_COUNT.map { sample, index, dedup_count -> [ [sample, index], dedup_count ] }
    .set { RNASEQ_INDIVIDUAL_GENOME_NODEDUP_COUNT_INDEXED }

    RNASEQ_CLIP_LOG_INDEXED_FOR_GENOME.join(RNASEQ_FILTER_LOG_INDEXED_FOR_GENOME)
            .join(RNASEQ_GENOME_ALIGNMENT_LOG_INDEXED)
            .join(RNASEQ_GENOME_QPASS_COUNTS_INDEXED)
            .join(RNASEQ_INDIVIDUAL_GENOME_NODEDUP_COUNT_INDEXED)
            .flatten()
            .collate(7)
            .set { RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }

    process rnaseq_individual_genome_alignment_stats {
        storeDir get_storedir('stats', true) + '/genome/individual'

        input:
        set val(sample), val(index), file(clip_log), file(filter_log), \
            file(genome_log), file(qpass_count), file(dedup_count) \
            from RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT

        output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome_individual.csv") \
            into RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS

        """
        rfc compile-step-stats \
            -n ${sample}.${index} \
            -c ${clip_log} \
            -f ${filter_log} \
            -g ${genome_log} \
            -q ${qpass_count} \
            -d ${dedup_count} \
            --label-prefix genome \
            -o ${sample}.${index}.rnaseq_genome_individual.csv
        """
    }

    // Group RNA-seq genome alignment outputs for merging (use merge channel)
    RNASEQ_GENOME_ALIGNMENT_BAM_FOR_MERGE.map { sample, index, bam -> [sample, bam] }.groupTuple()
    .set { RNASEQ_GENOME_ALIGNMENT_GROUPED_BAM }

    RNASEQ_GENOME_ALIGNMENT_ALIGNED.map { sample, index, fastq -> [sample, fastq] }.groupTuple()
    .set { RNASEQ_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ }

    RNASEQ_GENOME_ALIGNMENT_UNALIGNED.map { sample, index, fastq -> [sample, fastq] }.groupTuple()
    .set { RNASEQ_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ }

    RNASEQ_GENOME_ALIGNMENT_LOG_FOR_MERGE.map { sample, index, log -> [sample, log] }.groupTuple()
    .set { RNASEQ_GENOME_ALIGNMENT_GROUPED_LOG }

    // Join all grouped channels for merging
    RNASEQ_GENOME_ALIGNMENT_GROUPED_BAM.join(RNASEQ_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ)
                            .join(RNASEQ_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ)
                            .join(RNASEQ_GENOME_ALIGNMENT_GROUPED_LOG)
                            .set { RNASEQ_GENOME_ALIGNMENT_GROUPED_JOINT }

    process rnaseq_merge_genome_alignment {
        storeDir get_storedir('genome_alignment', true) + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(bam), file(aligned_fastq), \
          file(unaligned_fastq), file(alignment_log) from RNASEQ_GENOME_ALIGNMENT_GROUPED_JOINT

    output:
        set val(sample), file("${sample}.rnaseq_genome.bam") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_BAM
        set val(sample), file("${sample}.rnaseq_genome.bam.bai") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_BAI
        set val(sample), file("${sample}.rnaseq_genome.aligned.fastq.gz") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_ALIGNED_FASTQ
        set val(sample), file("${sample}.rnaseq_genome.unaligned.fastq.gz") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_UNALIGNED_FASTQ
        set val(sample), file("${sample}.rnaseq_genome.log") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_LOG
        set val(sample), file("${sample}.rnaseq_genome.csv") \
                      into RNASEQ_GENOME_ALIGNMENT_MERGED_CSV

        """
    samtools merge ${sample}.rnaseq_genome.bam ${bam} && samtools index ${sample}.rnaseq_genome.bam && \\
    zcat ${aligned_fastq} | gzip -c > ${sample}.rnaseq_genome.aligned.fastq.gz && \\
    zcat ${unaligned_fastq} | gzip -c > ${sample}.rnaseq_genome.unaligned.fastq.gz && \\
    cat ${alignment_log} > ${sample}.rnaseq_genome.log && \\
    rfc hisat2-log-to-csv \\
                       -l ${sample}.rnaseq_genome.log -n ${sample} -p rnaseq_genome \\
                       -o ${sample}.rnaseq_genome.csv
    """
    }

///////////////////////////////////////////////////////////////////////////////
//          R N A - S E Q   G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////

    // Split RNA-seq genome BAM channel for deduplication and bigwig generation
    RNASEQ_GENOME_ALIGNMENT_MERGED_BAM.into {
        RNASEQ_GENOME_BAM_FOR_DEDUP
        RNASEQ_GENOME_BAM_NODEDUP
        RNASEQ_GENOME_BAM_FOR_BIGWIG_REF
    }

    // RNA-seq genome deduplication now follows transcriptome pattern (BED-based instead of BAM-based)

    // RNA-seq genome deduplication (position-based only) - now uses merged BED approach
    process rnaseq_genome_deduplicate {
        storeDir get_storedir('deduplication', true) + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(bed) from RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP

        output:
        set val(sample), file("${sample}.rnaseq_genome.post_dedup.bed") \
            into RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP

        when:
        rnaseq_dedup_method == 'position'

        """
        rfc dedup -i ${bed} -o ${sample}.rnaseq_genome.post_dedup.bed
        """
    }

    // Split dedup BED channel for separation and BAM conversion
    RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.into {
        RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_SEP
        RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_BAM
    }

    // Individual count files will be created by the separation process

    // Create channel for genome separation (following transcriptome pattern)
    RNASEQ_GENOME_BED_FOR_INDEX_SEP_PRE
    .map { sample, index, file -> [sample, index] }
    .combine(RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_SEP, by:0)
    .set { RNASEQ_GENOME_BED_FOR_INDEX_SEP_POST_DEDUP }

    // Separate RNA-seq genome merged post_dedup bed back into individual lane files
    process rnaseq_separate_genome_bed_post_dedup {
        storeDir  get_storedir('deduplication', true) + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(bed) from RNASEQ_GENOME_BED_FOR_INDEX_SEP_POST_DEDUP

    output:
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.post_dedup.bed") \
     into RNASEQ_GENOME_BED_DEDUPLICATED
        set val(sample), val(index), file("${sample}.${index}.rnaseq_genome.count_after_dedup.txt")\
     into RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP_SEPARATED

    when:
    rnaseq_dedup_method == 'position'

        """
    awk -v this_sample=${sample}.${index} \
     '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} > ${sample}.${index}.rnaseq_genome.post_dedup.bed \
      && wc -l ${sample}.${index}.rnaseq_genome.post_dedup.bed > ${sample}.${index}.rnaseq_genome.count_after_dedup.txt
    """
    }

    // Select the appropriate count channel based on RNA-seq deduplication method
    if (rnaseq_dedup_method == 'position') {
        RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP_SEPARATED.set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT }
    }
    else {
        RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP.set { RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT }
    }

///////////////////////////////////////////////////////////////////////////////
//  E N D   R N A - S E Q   B I G W I G   G E N E R A T I O N  (INDIVIDUAL)
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  R N A - S E Q   B I G W I G   G E N E R A T I O N  (MERGED)
///////////////////////////////////////////////////////////////////////////////

    // Convert deduplicated BED to BAM for bigWig generation
    // Combine dedup BED with merged BAM (to get genome info)
    RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_BAM
    .combine(RNASEQ_GENOME_BAM_FOR_BIGWIG_REF, by:0)
    .set { RNASEQ_DEDUP_BED_WITH_BAM }

    process rnaseq_dedup_bed_to_bam {
        storeDir get_storedir('deduplication', true) + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(bed), file(ref_bam) from RNASEQ_DEDUP_BED_WITH_BAM

        output:
        set val(sample), file("${sample}.rnaseq_genome.post_dedup.bam"), \
                         file("${sample}.rnaseq_genome.post_dedup.bam.bai") \
            into RNASEQ_GENOME_DEDUP_BAM_FOR_BIGWIG

        when:
        rnaseq_dedup_method == 'position'

        """
        samtools view -H ${ref_bam} | grep '@SQ' | sed 's/@SQ\\tSN://g' | sed 's/\\tLN:/\\t/g' > genome.sizes
        samtools view -H ${ref_bam} > header.sam
        bedtools bedtobam -i ${bed} -g genome.sizes > unsorted.bam
        samtools sort -@ ${task.cpus} -o sorted.bam unsorted.bam
        samtools reheader header.sam sorted.bam > ${sample}.rnaseq_genome.post_dedup.bam
        samtools index -@ ${task.cpus} ${sample}.rnaseq_genome.post_dedup.bam
        """
    }

    // Generate strand-specific bigWig files from deduplicated BAMs
    process rnaseq_create_strand_specific_bigwigs {
        storeDir get_storedir('deduplication', true) + '/bigwigs/' + params.output.merged_lane_directory

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_bigwig'

        input:
        set val(sample), file(bam), file(bai) from RNASEQ_GENOME_DEDUP_BAM_FOR_BIGWIG

        output:
        set val(sample), file("${sample}.rnaseq.dedup.plus.bigWig"), \
                         file("${sample}.rnaseq.dedup.minus.bigWig") \
            into RNASEQ_STRAND_SPECIFIC_BIGWIGS

        when:
        rnaseq_dedup_method == 'position'

        """
        bamCoverage -b ${bam} -o ${sample}.rnaseq.dedup.plus.bigWig \
            --filterRNAstrand forward --binSize 1 -p ${task.cpus}

        bamCoverage -b ${bam} -o ${sample}.rnaseq.dedup.minus.bigWig \
            --filterRNAstrand reverse --binSize 1 -p ${task.cpus}
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

    input:
        file(stat_table) from RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

    output:
        file('rnaseq_genome_individual_essential.csv') \
          into COMBINED_INDIVIDUAL_RNASEQ_GENOME_ALIGNMENT_STATS

        script:
        if (stat_table.size() == 0) {
            // Create empty stats file when no input is available
            """
            echo "No RNA-seq individual statistics data available" > rnaseq_genome_individual_essential.csv
            echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> rnaseq_genome_individual_essential.csv
            """
        } else {
            """
      rfc merge overall-stats \
       -o raw_combined_individual_rnaseq_genome_aln_stats.csv \
          ${stat_table} && \
      rfc stats-percentage \
       -i raw_combined_individual_rnaseq_genome_aln_stats.csv \
       -l genome \
       -o rnaseq_genome_individual_essential.csv
        """
        }
    }

    // Sum individual RNA-seq genome stats by sample
    RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GROUPING
    .map { sample, index, file -> [ sample, file ] }
    .groupTuple()
    .set { RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED }

    process sum_individual_rnaseq_genome_alignment_stats {
        executor 'local'
        storeDir get_storedir('stats', true) + '/genome/merged'

        input:
        set val(sample), file(stat_files) from RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED

        output:
        set val(sample), file("${sample}.rnaseq_genome_merged.csv") into RNASEQ_GENOME_MERGED_ALIGNMENT_STATS

        """
        rfc sum-stats -n ${sample} -o ${sample}.rnaseq_genome_merged.csv ${stat_files}
        """
    }

    // Combine merged RNA-seq genome stats
    RNASEQ_GENOME_MERGED_ALIGNMENT_STATS
    .map { sample, stats_file -> stats_file }
    .toSortedList().set { RNASEQ_GENOME_MERGED_ALIGNMENT_STATS_COLLECTED }

    process combine_merged_rnaseq_genome_alignment_stats {
        executor 'local'
        storeDir get_storedir('stats', true) + '/genome/merged'

    input:
        file(stat_table) from RNASEQ_GENOME_MERGED_ALIGNMENT_STATS_COLLECTED

    output:
        file('rnaseq_genome_merged_essential.csv') \
          into COMBINED_MERGED_RNASEQ_GENOME_ALIGNMENT_STATS

        script:
        if (stat_table.size() == 0) {
            // Create empty stats file when no input is available
            """
            echo "No RNA-seq merged statistics data available" > rnaseq_genome_merged_essential.csv
            echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> rnaseq_genome_merged_essential.csv
            """
        } else {
            """
      rfc merge overall-stats \
       -o raw_combined_merged_rnaseq_genome_aln_stats.csv \
          ${stat_table} && \
      rfc stats-percentage \
       -i raw_combined_merged_rnaseq_genome_aln_stats.csv \
       -l genome \
       -o rnaseq_genome_merged_essential.csv
        """
        }
    }

///////////////////////////////////////////////////////////////////////////////
//          E N D   R N A - S E Q   G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////
} // if (do_rnaseq)
// RNA-Seq
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Merge Ribos*/

if (do_create_ribo) {
    if (do_rnaseq) {
        RIBO_WITH_RNASEQ.into { RIBO_FOR_MERGE_PRE; RIBO_FOR_COUNT }
    }
  else {
        RIBO_AFTER_CREATION.into { RIBO_FOR_MERGE_PRE; RIBO_FOR_COUNT }
  }
    RIBO_FOR_MERGE_PRE.map { sample, ribo -> [ribo] }.flatten().collect()
                    .set { RIBO_FOR_MERGE }

    process merge_ribos {
        publishDir get_publishdir('ribo'), mode:'copy'

    input:
        file(sample_ribo) from RIBO_FOR_MERGE
        val(ribo_count) from RIBO_FOR_COUNT.count()

    output:
        file('all.ribo') into ALL_RIBO

    script:
        if (ribo_count > 1) {
            command = "ribopy merge all.ribo ${sample_ribo}"
    } else {
            command = "ln -s ${sample_ribo} all.ribo"
        }

        """
    ${command}
    """
    }
}

// Merge Ribos
////////////////////////////////////////////////////////////////////////////////
