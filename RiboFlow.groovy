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
            ribo: 'ribo'
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
if (!params.containsKey('do_ribo_creation')) {
    params.do_ribo_creation = false
}
if (!params.containsKey('do_metadata')) {
    params.do_metadata = false
}

rnaseq_dedup_method = params.get('rnaseq', [:]).get('dedup_method', 'position').toString()
rnaseq_library_strandedness = params.get('rnaseq', [:]).get('library_strandedness', 'reverse').toString()

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
     ${params.umi_tools_extract_arguments ?: ''}
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
            | samtools view -b - \
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

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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
    set -o pipefail
    hisat2 ${params.alignment_arguments.genome} \
           -x ${genome_base} -U ${fastq} \
           -p ${task.cpus} \
           --no-softclip \
           --al-gz ${sample}.${index}.genome_alignment.aligned.fastq.gz \
           --un-gz ${sample}.${index}.genome_alignment.unaligned.fastq.gz \
               2> ${sample}.${index}.genome_alignment.log \
           | samtools view -b - \
           | samtools addreplacerg -r "ID:${sample}.${index}" -r "SM:${sample}" -r "PL:ILLUMINA" - \
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
/* GENOME ALIGNMENT FASTQC */

process genome_aligned_fastqc {
    publishDir get_publishdir("fastqc") + "/genome_aligned", mode: 'copy'

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
    
    # Check if file is empty (less than or equal to 20 bytes, which is just the gzip header)
    # Use -L to dereference symlink to get actual file size
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
    publishDir get_publishdir("fastqc") + "/genome_unaligned", mode: 'copy'

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
    
    # Check if file is empty (less than or equal to 20 bytes, which is just the gzip header)
    # Use -L to dereference symlink to get actual file size
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
    samtools view -bq ${params.mapping_quality_cutoff} ${bam} > ${sample}.${index}.genome_alignment.qpass.bam && \
    samtools index ${sample}.${index}.genome_alignment.qpass.bam
    samtools view -c ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome.qpass.count
    samtools flagstat ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome_alignment.qpass.stats
    """
}

GENOME_QPASS_COUNTS_OUT.into {
    GENOME_QPASS_COUNTS
    GENOME_QPASS_COUNTS_FOR_PSITE
    GENOME_QPASS_COUNTS_FOR_STATS
}

def psite_offset_file_exists = false
if (params.containsKey('psite_offset') && params.psite_offset.containsKey('offset_file')) {
    psite_offset_file_exists = file(params.psite_offset.offset_file).exists()
}

// Helper function to get GSM ID for a specific sample and lane index
// Supports both old experiment_mapping (backward compatibility) and new sample_matching structure
def get_gsm_id_for_lane(sample, index, file_exists_check = null) {
    def file_exists = (file_exists_check != null) ? file_exists_check : psite_offset_file_exists
    if (!file_exists || !params.containsKey('psite_offset')) {
        return null
    }
    
    def psite_config = params.psite_offset
    
    if (psite_config.containsKey('sample_matching') && psite_config.sample_matching.containsKey(sample)) {
        def sample_config = psite_config.sample_matching[sample]
        if (sample_config.containsKey('lanes')) {
            def lanes = sample_config.lanes
            if (lanes.containsKey(index.toString())) {
                return lanes[index.toString()]
            } else if (lanes.containsKey(index)) {
                return lanes[index]
            }
        }
    }
    
    if (psite_config.containsKey('experiment_mapping') && psite_config.experiment_mapping.containsKey(sample)) {
        return psite_config.experiment_mapping[sample]
    }
    
    return null
}

if (psite_offset_file_exists) {
    GENOME_ALIGNMENT_QPASS_BAM.into {
        GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE
        GENOME_ALIGNMENT_QPASS_BAM_FOR_STATS
        GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_REF
        GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_NONE_COMBINE
    }
    GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_NONE_COMBINE
        .combine(GENOME_ALIGNMENT_QPASS_BAI, by: [0, 1])
        .into {
            GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_NONE
            GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_NONE_FOR_POST_DEDUP
        }
} else {
    GENOME_ALIGNMENT_QPASS_BAM.into {
        GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE
        GENOME_ALIGNMENT_QPASS_BAM_FOR_STATS
    }
    // Create empty channel for missing variable when psite_offset_file_exists is false
    Channel.empty().set { GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_NONE }
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
    dedup_method == 'umi_tools' || dedup_method == 'none' || dedup_method == 'position'

        """
  if [ "${dedup_method}" == "umi_tools" ] || [ "${dedup_method}" == "position" ]; then
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

if (psite_offset_file_exists && dedup_method == 'position') {
    GENOME_MERGED_BAM_QPASS_PRE_DEDUP_RAW.into {
        GENOME_MERGED_BAM_QPASS_PRE_DEDUP_MAIN
        GENOME_MERGED_BAM_QPASS_FOR_PSITE
    }
} else {
    GENOME_MERGED_BAM_QPASS_PRE_DEDUP_RAW.set { GENOME_MERGED_BAM_QPASS_PRE_DEDUP_MAIN }
}

GENOME_MERGED_BAM_QPASS_PRE_DEDUP_MAIN.into {
    GENOME_MERGED_BAM_FOR_UMI_DEDUP
    GENOME_MERGED_BAM_FOR_POSITION_BAM_CONV
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

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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
    (for f in ${aligned_fastq}; do if [ -s "\$f" ] && gzip -t "\$f" 2>/dev/null; then zcat "\$f"; fi; done) | gzip -c > ${sample}.genome.aligned.fastq.gz && \
    (for f in ${unaligned_fastq}; do if [ -s "\$f" ] && gzip -t "\$f" 2>/dev/null; then zcat "\$f"; fi; done) | gzip -c > ${sample}.genome.unaligned.fastq.gz && \
    rfc merge-hisat2-logs -o ${sample}.genome.log ${alignment_log} && \
    rfc hisat2-log-to-csv -l ${sample}.genome.log -n ${sample} -p genome \
                      -o ${sample}.genome.csv
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
    storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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

if (psite_offset_file_exists && dedup_method == 'position') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_PSITE
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_BAM_EXTRACTION
    }
    
    process generate_post_dedup_counts_from_bed_position_with_psite {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(post_dedup_bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_STATS

        output:
        set val(sample), file("${sample}.count_after_dedup.txt") into GENOME_POSITION_DEDUP_COUNTS_WITH_PSITE

        when:
        psite_offset_file_exists && dedup_method == 'position'

        """
        wc -l < ${post_dedup_bed} > ${sample}.count_after_dedup.txt
        """
    }
    
    GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS
        .mix(GENOME_POSITION_DEDUP_COUNTS_WITH_PSITE)
        .set { GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS }
    
} else if (dedup_method == 'position') {
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
        .set { GENOME_BED_FOR_POSITION_SEPARATION_NO_PSITE }
    
    process separate_genome_bed_post_dedup_no_psite {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(bed) from GENOME_BED_FOR_POSITION_SEPARATION_NO_PSITE
        
        output:
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") \
                into GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION_NO_PSITE
        
        when:
        dedup_method == 'position' && !psite_offset_file_exists
        
        """
        # Extract reads for this specific lane using sample index column (column 7)
        awk -v this_sample=${sample}.${index} '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} > ${sample}.${index}.genome.post_dedup.bed
        """
    }
    
    process count_individual_separated_bed_position_no_psite {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(bed) from GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION_NO_PSITE
        
        output:
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION_NO_PSITE
        
        when:
        dedup_method == 'position' && !psite_offset_file_exists
        
        """
        wc -l < ${bed} > ${sample}.${index}.count_after_dedup.txt
        """
    }
} else if (dedup_method == 'umi_tools') {
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX }
} else {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_RAW.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX }
}

process genome_deduplicate_umi_tools {
    storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam), file(bai) from GENOME_MERGED_BAM_FOR_UMI_DEDUP

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
    storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

    input:
    set val(sample), file(bam) from GENOME_UMI_DEDUP_BAM_FOR_BED

    output:
    set val(sample), file("${sample}.genome.post_dedup.bed") into GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI

    when:
    dedup_method == 'umi_tools'

    """
    bamToBed -i ${bam} > ${sample}.genome.post_dedup.bed
    """
}

if (dedup_method == 'umi_tools') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_INDEX
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_STATS
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_MIX
    }
} else {
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_MIX }
}

if (dedup_method == 'umi_tools' && psite_offset_file_exists) {
    GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED
        .map { sample, bed -> [sample, bed] }
        .groupTuple()
        .set { GENOME_INDIVIDUAL_BED_WITH_INDEX_GROUPED }
    
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_INDEX
        .combine(GENOME_INDIVIDUAL_BED_WITH_INDEX_GROUPED, by: 0)
        .set { GENOME_BED_FOR_ADDING_SAMPLE_INDEX_UMI }
    
    process add_sample_index_to_merged_umi_dedup_bed {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory
        
        input:
            set val(sample), file(merged_bed), file(individual_beds) from GENOME_BED_FOR_ADDING_SAMPLE_INDEX_UMI
        
        output:
            set val(sample), file("${sample}.genome.post_dedup.with_sample_index.bed") \
                into GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_WITH_INDEX
        
        when:
        dedup_method == 'umi_tools' && psite_offset_file_exists
        
        """
        # Create a lookup table mapping read names (column 4) to sample.index (column 7)
        # from individual BED files (which have sample index in column 7)
        for bed_file in ${individual_beds}; do
            awk '{print \$4"\\t"\$7}' \$bed_file
        done | sort -u > read_to_sample_index.map
        
        # Match merged BED read names against lookup table and add sample index column
        awk 'BEGIN {
            while ((getline line < "read_to_sample_index.map") > 0) {
                split(line, parts, "\\t")
                read_name = parts[1]
                sample_index = parts[2]
                lookup[read_name] = sample_index
            }
            close("read_to_sample_index.map")
        }
        {
            read_name = \$4
            if (read_name in lookup) {
                print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"lookup[read_name]
            } else {
                # If read name not found, warn but still output (some reads may be removed during dedup)
                print "WARNING: Read name " read_name " not found in lookup" > "/dev/stderr"
                print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\tUNKNOWN"
            }
        }' ${merged_bed} > ${sample}.genome.post_dedup.with_sample_index.bed
        """
        }
}

GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_MIX.mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_FOR_MIX)
.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP }

if (dedup_method == 'umi_tools') {
    GENOME_UMI_TOOLS_DEDUP_BAM.into {
        GENOME_UMI_TOOLS_DEDUP_BAM_FOR_SPLITTING
        GENOME_UMI_TOOLS_DEDUP_BAM_FOR_OTHER
    }
}

///////////////////////////////////////////////////////////////////////////////
//          P - S I T E   O F F S E T   C O R R E C T I O N
///////////////////////////////////////////////////////////////////////////////

if (dedup_method == 'umi_tools') {
    GENOME_UMI_TOOLS_DEDUP_BAM_FOR_SPLITTING
        .cross(GENOME_INDIVIDUAL_BED_FOR_SPLITTING.map { sample, index, bed -> [sample, index] }.unique())
        .map { merged_data, individual_data ->
            [individual_data[0], individual_data[1], merged_data[1]]
        }
        .set { GENOME_DEDUP_BAM_FOR_SPLITTING }

    process split_genome_dedup_bam_to_individual {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(merged_bam) from GENOME_DEDUP_BAM_FOR_SPLITTING

        output:
            set val(sample), val(index), file("${sample}.${index}.post_dedup.bam"), \
                             file("${sample}.${index}.post_dedup.bam.bai") into GENOME_INDIVIDUAL_DEDUP_BAM
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING_UMI
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") into GENOME_INDIVIDUAL_POST_DEDUP_BED_UMI

        when:
        dedup_method == 'umi_tools'

        """
        # Check if read groups exist in the merged BAM
        if ! samtools view -H ${merged_bam} | grep -q "^@RG"; then
            echo "ERROR: No read groups found in merged BAM ${merged_bam}" >&2
            echo "Read groups are required for splitting merged UMI dedup BAM" >&2
            exit 1
        fi
        
        # Split by read group ID
        samtools view -B -r ${sample}.${index} ${merged_bam} -o ${sample}.${index}.post_dedup.bam
        
        # Check if splitting produced any reads
        read_count=\$(samtools view -c ${sample}.${index}.post_dedup.bam)
        
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
    
    GENOME_INDIVIDUAL_DEDUP_BAM.into {
        GENOME_INDIVIDUAL_DEDUP_BAM_FOR_PSITE
        GENOME_INDIVIDUAL_DEDUP_BAM_FOR_STATS
    }
} else {
    Channel.empty().set { GENOME_INDIVIDUAL_DEDUP_BAM_FOR_PSITE }
    Channel.empty().set { GENOME_INDIVIDUAL_DEDUP_BAM_FOR_STATS }
}

if (psite_offset_file_exists) {
    def offsetPath = file(params.psite_offset.offset_file)
    Channel.fromPath(offsetPath.toString()).first().set { PSITE_OFFSET_FILE_UMI }
    Channel.fromPath(offsetPath.toString()).first().set { PSITE_OFFSET_FILE_POSITION }
    Channel.fromPath(offsetPath.toString()).first().set { PSITE_OFFSET_FILE_NONE }
} else {
    Channel.empty().set { PSITE_OFFSET_FILE_UMI }
    Channel.empty().set { PSITE_OFFSET_FILE_POSITION }
    Channel.empty().set { PSITE_OFFSET_FILE_NONE }
}

if (psite_offset_file_exists && dedup_method == 'position') {
    GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_REF
    .map { sample, index, bam -> [sample, bam] }
    .groupTuple()
    .map { sample, bams -> [sample, bams[0]] }
    .set { GENOME_QPASS_BAM_SINGLE_FOR_PSITE }
}


if (psite_offset_file_exists && dedup_method == 'umi_tools') {
    
    GENOME_QPASS_COUNTS_FOR_PSITE
        .map { sample, index, count_file -> [sample, index] }
        .unique()
        .combine(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI_WITH_INDEX, by: 0)
        .set { GENOME_BED_FOR_UMI_PSITE_SEPARATION }
    
    process separate_genome_bed_post_dedup_for_psite_umi {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(bed) from GENOME_BED_FOR_UMI_PSITE_SEPARATION
        
        output:
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") \
                into GENOME_INDIVIDUAL_BED_POST_DEDUP_UMI
        
        when:
        psite_offset_file_exists && dedup_method == 'umi_tools'
        
        """
        # Extract reads for this specific lane using sample index column (column 7)
        awk -v this_sample=${sample}.${index} '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} > ${sample}.${index}.genome.post_dedup.bed
        """
    }
    
    GENOME_INDIVIDUAL_BED_POST_DEDUP_UMI.into {
        GENOME_INDIVIDUAL_BED_POST_DEDUP_UMI_FOR_PSITE
        GENOME_INDIVIDUAL_BED_POST_DEDUP_UMI_FOR_COUNT
    }
    
    process apply_psite_correction_umi_bed_individual {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(dedup_bed) from GENOME_INDIVIDUAL_BED_POST_DEDUP_UMI_FOR_PSITE
            file(offset_csv) from PSITE_OFFSET_FILE_UMI

        output:
            set val(sample), val(index), file("${sample}.${index}.psite.bed") \
                into GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_UMI

        when:
        psite_offset_file_exists && dedup_method == 'umi_tools'

        script:
        def experiment_id = get_gsm_id_for_lane(sample, index, true)

        if (experiment_id == null) {
            error "No sample_matching found for sample '${sample}' lane ${index} in psite_offset configuration"
        }

        """
        python3 ${workflow.projectDir}/scripts/apply_psite_offsets_bed.py \\
            -i ${dedup_bed} \\
            -o ${sample}.${index}.psite.bed \\
            -c ${offset_csv} \\
            -e ${experiment_id} \\
            -s ${sample}.${index}
        """
    }
    
    GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_UMI.into {
        GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_UMI_FOR_COUNT
        GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_UMI_FOR_GROUPING
    }
    
    process count_individual_psite_bed_umi {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(psite_bed) from GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_UMI_FOR_COUNT
        
        output:
            set val(sample), val(index), file("${sample}.${index}.psite_count.txt") into GENOME_INDIVIDUAL_PSITE_COUNTS_UMI
        
        when:
        dedup_method == 'umi_tools' && psite_offset_file_exists
        
        """
        wc -l < ${psite_bed} > ${sample}.${index}.psite_count.txt
        """
    }
    
    process count_individual_separated_bed_umi {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(bed) from GENOME_INDIVIDUAL_BED_POST_DEDUP_UMI_FOR_COUNT
        
        output:
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION_UMI
        
        when:
        dedup_method == 'umi_tools' && psite_offset_file_exists
        
        """
        wc -l < ${bed} > ${sample}.${index}.count_after_dedup.txt
        """
    }
}

if (psite_offset_file_exists && dedup_method == 'position') {
    GENOME_QPASS_COUNTS_FOR_PSITE
        .map { sample, index, count_file -> [sample, index] }
        .unique()
        .combine(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_PSITE, by: 0)
        .set { GENOME_BED_FOR_POSITION_PSITE_SEPARATION }
    
    process separate_genome_bed_post_dedup_for_psite {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(bed) from GENOME_BED_FOR_POSITION_PSITE_SEPARATION
        
        output:
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") \
                into GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION
        
        when:
        psite_offset_file_exists && dedup_method == 'position'
        
        """
        # Extract reads for this specific lane using sample index column (column 7)
        awk -v this_sample=${sample}.${index} '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} > ${sample}.${index}.genome.post_dedup.bed
        """
    }
    
    GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION.into {
        GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION_FOR_PSITE
        GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION_FOR_COUNT
    }
    
    process apply_psite_correction_position_bed_individual {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(dedup_bed) from GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION_FOR_PSITE
            file(offset_csv) from PSITE_OFFSET_FILE_POSITION

        output:
            set val(sample), val(index), file("${sample}.${index}.psite.bed") \
                into GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_POSITION

        when:
        psite_offset_file_exists && dedup_method == 'position'

        script:
        def experiment_id = get_gsm_id_for_lane(sample, index, true)

        if (experiment_id == null) {
            error "No sample_matching found for sample '${sample}' lane ${index} in psite_offset configuration"
        }

        """
        python3 ${workflow.projectDir}/scripts/apply_psite_offsets_bed.py \\
            -i ${dedup_bed} \\
            -o ${sample}.${index}.psite.bed \\
            -c ${offset_csv} \\
            -e ${experiment_id} \\
            -s ${sample}.${index}
        """
    }
    
    process count_individual_separated_bed_position {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(bed) from GENOME_INDIVIDUAL_BED_POST_DEDUP_POSITION_FOR_COUNT
        
        output:
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION
        
        when:
        psite_offset_file_exists && dedup_method == 'position'
        
        """
        wc -l < ${bed} > ${sample}.${index}.count_after_dedup.txt
        """
    }
    
    GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_POSITION.into {
        GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_POSITION_FOR_COUNT
        GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_POSITION_FOR_GROUPING
    }
    
    process count_individual_psite_bed_position {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
        
        input:
            set val(sample), val(index), file(psite_bed) from GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_POSITION_FOR_COUNT
        
        output:
            set val(sample), val(index), file("${sample}.${index}.psite_count.txt") into GENOME_INDIVIDUAL_PSITE_COUNTS_POSITION
        
        when:
        dedup_method == 'position' && psite_offset_file_exists
        
        """
        wc -l < ${psite_bed} > ${sample}.${index}.psite_count.txt
        """
    }
    
    GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_POSITION_FOR_GROUPING
        .map { sample, index, bed -> [sample, bed] }
        .groupTuple()
        .set { GENOME_PSITE_CORRECTED_BED_POSITION_GROUPED }
}

process apply_psite_correction_umi_individual {
    storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory

    input:
        set val(sample), val(index), file(dedup_bam), file(dedup_bai) from GENOME_INDIVIDUAL_DEDUP_BAM_FOR_PSITE
        file(offset_csv) from PSITE_OFFSET_FILE_UMI

    output:
        set val(sample), val(index), file("${sample}.${index}.psite.bam"), file("${sample}.${index}.psite.bam.bai") into GENOME_INDIVIDUAL_PSITE_CORRECTED_BAM_UMI

    when:
    psite_offset_file_exists && dedup_method == 'umi_tools'

    script:
    def experiment_id = get_gsm_id_for_lane(sample, index, true)

    if (experiment_id == null) {
        error "No sample_matching found for sample '${sample}' lane ${index} in psite_offset configuration"
    }

    """
    python3 ${workflow.projectDir}/scripts/apply_psite_offsets.py \\
        -i ${dedup_bam} \\
        -o ${sample}.${index}.psite.bam \\
        -c ${offset_csv} \\
        -e ${experiment_id} \\
        -s ${sample}.${index} \\
        --index
    """
}

GENOME_INDIVIDUAL_PSITE_CORRECTED_BAM_UMI
    .map { sample, index, bam, bai -> [sample, bam] }
    .groupTuple()
    .set { GENOME_PSITE_CORRECTED_BAM_UMI_GROUPED }

Channel.empty().set { GENOME_PSITE_CORRECTED_BED_UMI_GROUPED }

if (dedup_method == 'umi_tools' && psite_offset_file_exists) {
    GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_UMI_FOR_GROUPING
        .map { sample, index, bed -> [sample, bed] }
        .groupTuple()
        .set { GENOME_PSITE_CORRECTED_BED_UMI_GROUPED }
}

process apply_psite_correction_none_individual {
    storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
    
    input:
        set val(sample), val(index), file(qpass_bam), file(qpass_bai) from GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_NONE
        file(offset_csv) from PSITE_OFFSET_FILE_NONE
    
    output:
        set val(sample), val(index), file("${sample}.${index}.psite.bam"), file("${sample}.${index}.psite.bam.bai") \
            into GENOME_INDIVIDUAL_PSITE_CORRECTED_BAM_NONE
        set val(sample), val(index), file("${sample}.${index}.psite.bed") \
            into GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_NONE
    
    when:
    psite_offset_file_exists && dedup_method == 'none'
    
    script:
    def experiment_id = get_gsm_id_for_lane(sample, index, true)
    if (experiment_id == null) {
        error "No sample_matching found for sample '${sample}' lane ${index} in psite_offset configuration"
    }
    """
    python3 ${workflow.projectDir}/scripts/apply_psite_offsets.py \\
        -i ${qpass_bam} \\
        -o ${sample}.${index}.psite.bam \\
        -c ${offset_csv} \\
        -e ${experiment_id} \\
        -s ${sample}.${index} \\
        --index

    bamToBed -i ${sample}.${index}.psite.bam > ${sample}.${index}.psite.bed
    """
}

GENOME_INDIVIDUAL_PSITE_CORRECTED_BAM_NONE
    .map { sample, index, bam, bai -> [sample, bam] }
    .groupTuple()
    .set { GENOME_PSITE_CORRECTED_BAM_NONE_GROUPED }

if (psite_offset_file_exists) {
    if (dedup_method == 'umi_tools') {
        GENOME_PSITE_CORRECTED_BAM_UMI_GROUPED.mix(GENOME_PSITE_CORRECTED_BAM_NONE_GROUPED)
            .set { GENOME_PSITE_CORRECTED_BAM_ALL_GROUPED }
    } else if (dedup_method == 'none') {
        GENOME_PSITE_CORRECTED_BAM_NONE_GROUPED.set { GENOME_PSITE_CORRECTED_BAM_ALL_GROUPED }
    } else {
        Channel.empty().set { GENOME_PSITE_CORRECTED_BAM_ALL_GROUPED }
    }
} else {
    Channel.empty().set { GENOME_PSITE_CORRECTED_BAM_ALL_GROUPED }
}

process merge_psite_corrected_bam {
    storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(bam_files) from GENOME_PSITE_CORRECTED_BAM_ALL_GROUPED

    output:
        set val(sample), file("${sample}.psite.bam"), file("${sample}.psite.bam.bai") into GENOME_PSITE_CORRECTED_BAM_MERGED

    when:
    psite_offset_file_exists && (dedup_method == 'umi_tools' || dedup_method == 'none')

    """
    samtools merge ${sample}.psite.bam ${bam_files}
    samtools index -@ ${task.cpus} ${sample}.psite.bam
    """
}

GENOME_PSITE_CORRECTED_BAM_MERGED.into {
    GENOME_PSITE_CORRECTED_BAM_MERGED_UMI
    GENOME_PSITE_CORRECTED_BAM_MERGED_NONE
}

if (psite_offset_file_exists && dedup_method == 'umi_tools') {
    GENOME_PSITE_CORRECTED_BAM_MERGED_UMI.set { GENOME_PSITE_CORRECTED_BAM_UMI_MERGED }
    GENOME_PSITE_CORRECTED_BAM_UMI_MERGED.set { GENOME_BAM_FOR_BIGWIG_FINAL }
} else if (psite_offset_file_exists && dedup_method == 'none') {
    GENOME_PSITE_CORRECTED_BAM_MERGED_NONE.set { GENOME_PSITE_CORRECTED_BAM_NONE_MERGED }
}

GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_NONE.into {
    GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_NONE_FOR_COUNT
    GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_NONE_FOR_GROUPING
}

process count_individual_psite_bed_none {
    storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory
    
    input:
        set val(sample), val(index), file(psite_bed) from GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_NONE_FOR_COUNT
    
    output:
        set val(sample), val(index), file("${sample}.${index}.psite_count.txt") into GENOME_INDIVIDUAL_PSITE_COUNTS_NONE
    
    when:
    dedup_method == 'none' && psite_offset_file_exists
    
    """
    wc -l < ${psite_bed} > ${sample}.${index}.psite_count.txt
    """
}

GENOME_INDIVIDUAL_PSITE_CORRECTED_BED_NONE_FOR_GROUPING
    .map { sample, index, bed -> [sample, bed] }
    .groupTuple()
    .set { GENOME_PSITE_CORRECTED_BED_NONE_GROUPED }

if (psite_offset_file_exists) {
    if (dedup_method == 'position') {
        GENOME_PSITE_CORRECTED_BED_POSITION_GROUPED
            .mix(GENOME_PSITE_CORRECTED_BED_UMI_GROUPED, GENOME_PSITE_CORRECTED_BED_NONE_GROUPED)
            .set { GENOME_PSITE_CORRECTED_BED_ALL_GROUPED }
    } else if (dedup_method == 'umi_tools') {
        GENOME_PSITE_CORRECTED_BED_UMI_GROUPED
            .mix(GENOME_PSITE_CORRECTED_BED_NONE_GROUPED)
            .set { GENOME_PSITE_CORRECTED_BED_ALL_GROUPED }
    } else if (dedup_method == 'none') {
        GENOME_PSITE_CORRECTED_BED_NONE_GROUPED.set { GENOME_PSITE_CORRECTED_BED_ALL_GROUPED }
    } else {
        Channel.empty().set { GENOME_PSITE_CORRECTED_BED_ALL_GROUPED }
    }
} else {
    Channel.empty().set { GENOME_PSITE_CORRECTED_BED_ALL_GROUPED }
}

process merge_psite_corrected_bed {
    storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

    input:
        set val(sample), file(bed_files) from GENOME_PSITE_CORRECTED_BED_ALL_GROUPED

    output:
        set val(sample), file("${sample}.psite.bed") into GENOME_PSITE_CORRECTED_BED_MERGED

    when:
    psite_offset_file_exists

    """
    cat ${bed_files} > ${sample}.psite.bed
    """
}

GENOME_PSITE_CORRECTED_BED_MERGED.into {
    GENOME_PSITE_CORRECTED_BED_MERGED_POSITION
    GENOME_PSITE_CORRECTED_BED_MERGED_UMI
    GENOME_PSITE_CORRECTED_BED_MERGED_NONE
}

if (psite_offset_file_exists && dedup_method == 'position') {
    GENOME_PSITE_CORRECTED_BED_MERGED_POSITION.into {
        GENOME_PSITE_CORRECTED_BED_MERGED_POSITION_FOR_STATS
    }
    GENOME_PSITE_CORRECTED_BED_MERGED_POSITION_FOR_STATS.set { GENOME_PSITE_BED_FOR_STATS }
    
    // Use the full-length deduplicated BED (before P-site correction) for BAM extraction
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_BAM_EXTRACTION
        .join(GENOME_MERGED_BAM_QPASS_FOR_PSITE)
        .map { sample, dedup_bed, qpass_bam, qpass_bai -> [sample, dedup_bed, qpass_bam, qpass_bai] }
        .set { GENOME_BED_FOR_DEDUP_WITH_QPASS_BAM_POSITION }
    
    process genome_extract_psite_reads_from_bed {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(dedup_bed), file(qpass_bam), file(qpass_bai) from GENOME_BED_FOR_DEDUP_WITH_QPASS_BAM_POSITION
        file(offset_csv) from PSITE_OFFSET_FILE_POSITION

        output:
        set val(sample), file("${sample}.psite.bam"), file("${sample}.psite.bam.bai") \
            into GENOME_PSITE_CORRECTED_BAM_POSITION_MERGED

        when:
        psite_offset_file_exists && dedup_method == 'position'
        
        script:
        def experiment_id = get_gsm_id_for_lane(sample, 'merged', true)
        
        if (experiment_id == null) {
             // For merged samples, try to get experiment ID from first constituent lane
             experiment_id = get_gsm_id_for_lane(sample, '1', true)
        }
        
        if (experiment_id == null) {
             // Try lane 1 as integer
             experiment_id = get_gsm_id_for_lane(sample, 1, true)
        }
        
        if (experiment_id == null) {
             // Try lane 2 as fallback
             experiment_id = get_gsm_id_for_lane(sample, '2', true)
        }
        
        if (experiment_id == null) {
             // Fallback: try looking up the sample directly in experiment_mapping
             if (params.psite_offset.experiment_mapping != null && params.psite_offset.experiment_mapping.containsKey(sample)) {
                 experiment_id = params.psite_offset.experiment_mapping[sample]
             }
        }
        
        if (experiment_id == null) {
            error "No experiment ID found for sample '${sample}' in psite_offset configuration. Tried: merged, '1', 1, '2'"
        }

        """
        python3 ${workflow.projectDir}/scripts/extract_reads_from_dedup_bed.py \\
            --bam ${qpass_bam} \\
            --bed ${dedup_bed} \\
            --output ${sample}.post_dedup.bam
            
        samtools index ${sample}.post_dedup.bam
        
        python3 ${workflow.projectDir}/scripts/apply_psite_offsets.py \\
            -i ${sample}.post_dedup.bam \\
            -o ${sample}.psite.bam \\
            -c ${offset_csv} \\
            -e ${experiment_id} \\
            -s ${sample} \\
            --index
        """
    }
    
    GENOME_PSITE_CORRECTED_BAM_POSITION_MERGED.set { GENOME_BAM_FOR_BIGWIG_FINAL }
    
    process generate_psite_counts_from_bed_position {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(psite_bed) from GENOME_PSITE_BED_FOR_STATS

        output:
        set val(sample), file("${sample}.psite_count.txt") into GENOME_PSITE_COUNTS_FROM_BED_POSITION_OUTPUT

        """
        wc -l < ${psite_bed} > ${sample}.psite_count.txt
        """
    }
}

if (psite_offset_file_exists && dedup_method == 'umi_tools') {
    GENOME_PSITE_CORRECTED_BED_MERGED_UMI.set { GENOME_PSITE_CORRECTED_BED_UMI_MERGED }
}

if (psite_offset_file_exists && dedup_method == 'none') {
    GENOME_PSITE_CORRECTED_BED_MERGED_NONE.set { GENOME_PSITE_CORRECTED_BED_NONE_MERGED }
    GENOME_PSITE_CORRECTED_BAM_NONE_MERGED.set { GENOME_PSITE_CORRECTED_BAM_NONE }
    GENOME_PSITE_CORRECTED_BED_NONE_MERGED.set { GENOME_PSITE_CORRECTED_BED_NONE }
}


if (psite_offset_file_exists && dedup_method == 'umi_tools') {
    GENOME_PSITE_CORRECTED_BED_UMI_MERGED.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_RAW_UMI
        GENOME_PSITE_BED_FOR_STATS
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_RAW_UMI.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_PROCESSES
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_UMI
        GENOME_BED_FOR_RIBO
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_RAW_UMI.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL }

    process generate_psite_counts_from_bed_umi {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(psite_bed) from GENOME_PSITE_BED_FOR_STATS

        output:
            set val(sample), file("${sample}.psite_count.txt") into GENOME_PSITE_COUNTS_FROM_BED_UMI_OUTPUT

        """
        wc -l < ${psite_bed} > ${sample}.psite_count.txt
        """
    }

    process generate_merged_dedup_count_umi_with_psite {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(post_dedup_bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_PROCESSES

        output:
            set val(sample), file("${sample}.count_after_dedup.txt") into GENOME_MERGED_DEDUP_COUNT_UMI_WITH_PSITE

        when:
        dedup_method == 'umi_tools' && psite_offset_file_exists

        """
        wc -l < ${post_dedup_bed} > ${sample}.count_after_dedup.txt
        """
    }
    
    GENOME_MERGED_DEDUP_COUNT_UMI_WITH_PSITE.set { GENOME_MERGED_DEDUP_COUNT_UMI }
    
} else if (!psite_offset_file_exists && dedup_method == 'position') {
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION_FOR_BAM
        .join(GENOME_MERGED_BAM_FOR_POSITION_BAM_CONV)
        .map { sample, dedup_bed, bam, bai -> [sample, dedup_bed, bam] }
        .set { GENOME_BED_FOR_BAM_CONVERSION_POSITION }

    process genome_convert_dedup_bed_to_bam_position {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
        set val(sample), file(dedup_bed), file(qpass_bam) from GENOME_BED_FOR_BAM_CONVERSION_POSITION

        output:
        set val(sample), file("${sample}.post_dedup.bam"), file("${sample}.post_dedup.bam.bai") \
            into GENOME_POST_DEDUP_BAM_POSITION_MERGED

        when:
        !psite_offset_file_exists && dedup_method == 'position'

        """
        python3 ${workflow.projectDir}/scripts/extract_reads_from_dedup_bed.py \
            --bam ${qpass_bam} \
            --bed ${dedup_bed} \
            --output ${sample}.post_dedup.bam && \
        samtools index ${sample}.post_dedup.bam
        """
    }

    GENOME_POST_DEDUP_BAM_POSITION_MERGED.set { GENOME_BAM_FOR_BIGWIG_FINAL }

} else if (psite_offset_file_exists && dedup_method == 'none') {
    GENOME_PSITE_CORRECTED_BAM_NONE.into {
        GENOME_BAM_FOR_BIGWIG_FINAL
        GENOME_PSITE_COUNTS_FROM_BAM
    }
    GENOME_PSITE_CORRECTED_BED_NONE.into {
        GENOME_PSITE_BED_FOR_STATS
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_RAW_NONE
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_RAW_NONE.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_PROCESSES
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_NONE
        GENOME_BED_FOR_RIBO
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_RAW_NONE.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL }
    
    process create_individual_post_dedup_bed_none {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.individual_lane_directory

        input:
            set val(sample), val(index), file(qpass_bam) from GENOME_ALIGNMENT_QPASS_BAM_FOR_PSITE_NONE_FOR_POST_DEDUP.map { sample, index, bam, bai -> [sample, index, bam] }

        output:
            set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bed") into GENOME_INDIVIDUAL_POST_DEDUP_BED_NONE
            set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_NONE

        when:
        psite_offset_file_exists && dedup_method == 'none'

        """
        if [ `samtools view -c ${qpass_bam}` -eq 0 ]; then
            touch ${sample}.${index}.genome.post_dedup.bed
            echo "0" > ${sample}.${index}.count_after_dedup.txt
        else
            bamToBed -i ${qpass_bam} > ${sample}.${index}.genome.post_dedup.bed
            wc -l < ${sample}.${index}.genome.post_dedup.bed > ${sample}.${index}.count_after_dedup.txt
        fi
        """
    }
    
    process create_merged_post_dedup_bed_none {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(qpass_bam), file(qpass_bai) from GENOME_MERGED_BAM_QPASS_NONE

        output:
            set val(sample), file("${sample}.genome.post_dedup.bed") into GENOME_MERGED_POST_DEDUP_BED_NONE
            set val(sample), file("${sample}.count_after_dedup.txt") into GENOME_MERGED_DEDUP_COUNT_NONE

        when:
        psite_offset_file_exists && dedup_method == 'none'

        """
        if [ `samtools view -c ${qpass_bam}` -eq 0 ]; then
            touch ${sample}.genome.post_dedup.bed
            echo "0" > ${sample}.count_after_dedup.txt
        else
            bamToBed -i ${qpass_bam} > ${sample}.genome.post_dedup.bed
            wc -l < ${sample}.genome.post_dedup.bed > ${sample}.count_after_dedup.txt
        fi
        """
    }

    process generate_psite_counts_from_bam_none {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(psite_bam), file(psite_bai) from GENOME_PSITE_COUNTS_FROM_BAM

        output:
            set val(sample), file("${sample}.psite_count.txt") into GENOME_PSITE_COUNTS_FROM_BAM_OUTPUT_NONE

        """
        samtools view -c ${psite_bam} > ${sample}.psite_count.txt
        """
    }

} else if (dedup_method == 'umi_tools') {
    GENOME_UMI_TOOLS_DEDUP_BAM_FOR_OTHER.into {
        GENOME_BAM_FOR_DOWNSTREAM
        GENOME_DEDUP_BAM_FOR_INDEXING
    }
    
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_FINAL
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_FINAL.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_PROCESSES
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_UMI_NO_PSITE
        GENOME_BED_FOR_RIBO
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_FINAL.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL }
    
    process generate_merged_dedup_count_umi {
        storeDir get_storedir('alignment_ribo') + '/' + params.output.merged_lane_directory

        input:
            set val(sample), file(post_dedup_bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT

        output:
            set val(sample), file("${sample}.count_after_dedup.txt") into GENOME_MERGED_DEDUP_COUNT_UMI

        when:
        dedup_method == 'umi_tools' && !psite_offset_file_exists

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
        !psite_offset_file_exists && dedup_method == 'umi_tools'

            """
        samtools index -@ ${task.cpus} ${bam}
        """
    }
} else {
    if (psite_offset_file_exists && dedup_method == 'none') {
    } else {
        if (!psite_offset_file_exists && dedup_method == 'none') {
            GENOME_MERGED_BAM_QPASS_NONE.set { GENOME_BAM_FOR_BIGWIG_FINAL }
        }
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.into {
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_PROCESSES
        GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_NODEDUP
        GENOME_BED_FOR_RIBO
    }
    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP.set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL }
    
    // No metadata channel needed for this branch as it's not using dedup or has psite
    Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_NODEDUP }
    Channel.empty().set { GENOME_BED_FOR_RIBO }
    }
}


if (dedup_method == 'umi_tools') {
    if (psite_offset_file_exists) {
        GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION_UMI.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    } else {
        GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING_UMI.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    }
}
else if (dedup_method == 'position') {
    if (psite_offset_file_exists) {
        GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    } else {
        GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SEPARATION_NO_PSITE.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    }
}
else {
    if (psite_offset_file_exists && dedup_method == 'none') {
        GENOME_INDIVIDUAL_DEDUP_COUNT_NONE.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    } else {
        GENOME_INDIVIDUAL_NODEDUP_COUNT.set { GENOME_INDIVIDUAL_DEDUP_COUNT }
    }
}

///////////////////////////////////////////////////////////////////////////////
//  R I B O - S E Q   B I G W I G   G E N E R A T I O N
///////////////////////////////////////////////////////////////////////////////

process genome_create_strand_specific_bigwigs {
    storeDir get_storedir('alignment_ribo') + '/bigwigs/' + params.output.merged_lane_directory

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_bigwig'

    input:
    set val(sample), file(bam), file(bai) from GENOME_BAM_FOR_BIGWIG_FINAL

    output:
    set val(sample), file("${sample}.psite.plus.bigWig"), \
                     file("${sample}.psite.minus.bigWig") \
        into GENOME_STRAND_SPECIFIC_BIGWIGS

    when:
    dedup_method == 'umi_tools' || dedup_method == 'none' || dedup_method == 'position'

    """
    if [ "${dedup_method}" == "none" ]; then
        # For nodedup ribo-seq, swap stranded information
        bamCoverage -b ${bam} -o ${sample}.psite.plus.bigWig \
            --filterRNAstrand reverse --binSize 1 -p ${task.cpus}
        bamCoverage -b ${bam} -o ${sample}.psite.minus.bigWig \
            --filterRNAstrand forward --binSize 1 -p ${task.cpus}
    else
        # For umi_tools, use standard stranded information
        bamCoverage -b ${bam} -o ${sample}.psite.plus.bigWig \
            --filterRNAstrand forward --binSize 1 -p ${task.cpus}
        bamCoverage -b ${bam} -o ${sample}.psite.minus.bigWig \
            --filterRNAstrand reverse --binSize 1 -p ${task.cpus}
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

if (psite_offset_file_exists) {
    if (dedup_method == 'umi_tools') {
        GENOME_PSITE_COUNTS_FROM_BED_UMI_OUTPUT
            .map { sample, psite_count_file -> [sample, psite_count_file] }
            .set { GENOME_PSITE_COUNTS_BY_SAMPLE_FOR_STATS }
    } else if (dedup_method == 'position') {
        GENOME_PSITE_COUNTS_FROM_BED_POSITION_OUTPUT
            .map { sample, psite_count_file -> [sample, psite_count_file] }
            .set { GENOME_PSITE_COUNTS_BY_SAMPLE_FOR_STATS }
    } else if (dedup_method == 'none') {
        GENOME_PSITE_COUNTS_FROM_BAM_OUTPUT_NONE
            .map { sample, psite_count_file -> [sample, psite_count_file] }
            .set { GENOME_PSITE_COUNTS_BY_SAMPLE_FOR_STATS }
    }
    
    if (dedup_method == 'umi_tools') {
        GENOME_INDIVIDUAL_PSITE_COUNTS_UMI
            .map { sample, index, count_file -> [[sample, index], count_file] }
            .set { GENOME_PSITE_COUNTS_INDEXED }
    } else if (dedup_method == 'position') {
        GENOME_INDIVIDUAL_PSITE_COUNTS_POSITION
            .map { sample, index, count_file -> [[sample, index], count_file] }
            .set { GENOME_PSITE_COUNTS_INDEXED }
    } else if (dedup_method == 'none') {
        GENOME_INDIVIDUAL_PSITE_COUNTS_NONE
            .map { sample, index, count_file -> [[sample, index], count_file] }
            .set { GENOME_PSITE_COUNTS_INDEXED }
    } else {
        GENOME_PSITE_COUNTS_BY_SAMPLE_FOR_STATS
            .join(GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED.map { key, dedup_file -> [key[0], [key, dedup_file]] }, by: 0)
            .map { sample, psite_count_file, key_and_dedup -> 
                def key = key_and_dedup[0]  // [sample, index] from original structure
                [key, psite_count_file] 
            }
            .set { GENOME_PSITE_COUNTS_INDEXED }
    }
}

CLIP_LOG_INDEXED_FOR_GENOME.join(FILTER_LOG_INDEXED_FOR_GENOME)
            .join(GENOME_ALIGNMENT_LOG_STATS_INDEXED)
            .join(GENOME_QPASS_COUNTS_INDEXED)
            .join(GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED)
            .map { key, clip_log, filter_log, genome_log, qpass_count, dedup_count -> 
                [key, clip_log, filter_log, genome_log, qpass_count, dedup_count] 
            }
            .set { GENOME_STATS_BASE_INPUT }

if (psite_offset_file_exists) {
    GENOME_STATS_BASE_INPUT
        .join(GENOME_PSITE_COUNTS_INDEXED, by: 0)
        .map { key, clip_log, filter_log, genome_log, qpass_count, dedup_count, psite_count_file -> 
            [key[0], key[1], clip_log, filter_log, genome_log, qpass_count, dedup_count, psite_count_file]
        }
        .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }
} else {
    GENOME_STATS_BASE_INPUT
        .map { key, clip_log, filter_log, genome_log, qpass_count, dedup_count -> 
            [key[0], key[1], clip_log, filter_log, genome_log, qpass_count, dedup_count, null]
        }
        .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }
}

process individual_genome_alignment_stats {
    storeDir get_storedir('stats') + '/genome/individual'

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

input:
    set val(sample), val(index), file(clip_log), file(filter_log),\
    file(genome_log), file(qpass_count),\
    file(dedup_count), val(psite_count_file)\
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
      -t ${genome_log} \\
      -q ${qpass_count} \\
      -d ${dedup_count} \\
      -l genome \\
      -o ${sample}.${index}.genome_individual.csv"
    # Add P-site count if available (passed through channel)
    if [[ -n "${params.psite_offset}" ]]; then
        if [ -n "${psite_count_file}" ] && [ "${psite_count_file}" != "null" ] && [ -f "${psite_count_file}" ]; then
            cmd="\${cmd} -p ${psite_count_file}"
        fi
    fi

    eval "\$cmd"
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

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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

GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GROUPING
  .map { sample, index, file -> [ sample, file ] }
  .groupTuple()
  .set { GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED }

if (dedup_method == 'position') {
    GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
        .combine(GENOME_POSITION_DEDUP_COUNTS_FOR_MERGED_STATS, by: 0)
        .set { GENOME_STATS_INPUT_WITH_DEDUP }
    if (psite_offset_file_exists) {
        GENOME_STATS_INPUT_WITH_DEDUP
            .combine(GENOME_PSITE_COUNTS_BY_SAMPLE_FOR_STATS, by: 0)
            .map { sample, stat_files, merged_dedup_count_file, merged_psite_count_file -> 
                [sample, stat_files, merged_dedup_count_file, merged_psite_count_file]
            }
            .set { GENOME_STATS_INPUT }
    } else {
        GENOME_STATS_INPUT_WITH_DEDUP
            .map { sample, stat_files, merged_dedup_count_file -> 
                [sample, stat_files, merged_dedup_count_file, null]
            }
            .set { GENOME_STATS_INPUT }
    }
} else if (dedup_method == 'umi_tools') {
    GENOME_MERGED_DEDUP_COUNT_UMI
        .map { sample, count_file -> [sample, count_file] }
        .set { GENOME_UMI_DEDUP_COUNTS_FOR_MERGED_STATS }
    GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
        .combine(GENOME_UMI_DEDUP_COUNTS_FOR_MERGED_STATS, by: 0)
        .set { GENOME_STATS_INPUT_WITH_DEDUP }
    if (psite_offset_file_exists) {
        GENOME_STATS_INPUT_WITH_DEDUP
            .combine(GENOME_PSITE_COUNTS_BY_SAMPLE_FOR_STATS, by: 0)
            .map { sample, stat_files, merged_dedup_count_file, merged_psite_count_file -> 
                [sample, stat_files, merged_dedup_count_file, merged_psite_count_file]
            }
            .set { GENOME_STATS_INPUT }
    } else {
        GENOME_STATS_INPUT_WITH_DEDUP
            .map { sample, stat_files, merged_dedup_count_file -> 
                [sample, stat_files, merged_dedup_count_file, null]
            }
            .set { GENOME_STATS_INPUT }
    }
} else {
    if (psite_offset_file_exists) {
        GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
            .combine(GENOME_PSITE_COUNTS_BY_SAMPLE_FOR_STATS, by: 0)
            .map { sample, stat_files, merged_psite_count_file -> 
                [sample, stat_files, null, merged_psite_count_file]
            }
            .set { GENOME_STATS_INPUT }
    } else {
        GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED
            .map { sample, stat_files -> 
                [sample, stat_files, null, null]
            }
            .set { GENOME_STATS_INPUT }
    }
}

process sum_individual_genome_alignment_stats {
    executor 'local'
    storeDir get_storedir('stats') + '/genome/merged'

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

input:
    // Use val() for count files - they come from process outputs and Nextflow will resolve their paths correctly
    // The script handles null/empty values by checking if the path string is 'null' or empty
    set val(sample), file(stat_files), val(merged_dedup_count_file), val(merged_psite_count_file) from GENOME_STATS_INPUT

output:
    set val(sample), file("${sample}.genome_merged.csv") into GENOME_MERGED_ALIGNMENT_STATS

    script:
    // Check if merged_dedup_count_file exists and is not empty
    merged_count_file_exists = merged_dedup_count_file != null && merged_dedup_count_file.toString() != 'null' && merged_dedup_count_file.toString() != ''
    // Check if merged_psite_count_file exists and is not empty
    merged_psite_count_file_exists = merged_psite_count_file != null && merged_psite_count_file.toString() != 'null' && merged_psite_count_file.toString() != ''
    
    if ((dedup_method == 'position' || dedup_method == 'umi_tools') && merged_count_file_exists) {
        if (merged_psite_count_file_exists) {
            """
            # Sum individual stats first (for all other metrics)
            rfc sum-stats -n ${sample} -o ${sample}.genome_merged.tmp.csv ${stat_files}
            
            # Update the merged dedup count and merged P-site count from the actual merged counts
            python3 ${workflow.projectDir}/scripts/update_merged_stats_with_counts.py \\
                --dedup-count-file ${merged_dedup_count_file} \\
                --psite-count-file ${merged_psite_count_file} \\
                --input-csv ${sample}.genome_merged.tmp.csv \\
                --output-csv ${sample}.genome_merged.csv
            """
        } else {
            """
            # Sum individual stats first (for all other metrics)
            rfc sum-stats -n ${sample} -o ${sample}.genome_merged.tmp.csv ${stat_files}
            
            # Update the merged dedup count from the actual merged count
            python3 ${workflow.projectDir}/scripts/update_merged_stats_with_counts.py \\
                --dedup-count-file ${merged_dedup_count_file} \\
                --input-csv ${sample}.genome_merged.tmp.csv \\
                --output-csv ${sample}.genome_merged.csv
            """
        }
    } else if (merged_psite_count_file_exists) {
        """
        # Sum individual stats first (for all other metrics)
        rfc sum-stats -n ${sample} -o ${sample}.genome_merged.tmp.csv ${stat_files}
        
        # Update the merged P-site count from the actual merged count
        python3 ${workflow.projectDir}/scripts/update_merged_stats_with_psite_only.py \\
            --psite-count-file ${merged_psite_count_file} \\
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

    beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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

// Initialize all metadata channels as empty (will be populated conditionally)
Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_UMI }
Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_POSITION }
Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_NONE }
Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_UMI_NO_PSITE }
Channel.empty().set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_NODEDUP }


// Combine all metadata channels
GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_UMI
    .mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_POSITION)
    .mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_NONE)
    .mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_UMI_NO_PSITE)
    .mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA_NODEDUP)
    .set { GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA }

if (do_metadata) {
    meta_base = params.input.metadata.get('base', '')
    if (meta_base != '' && !meta_base.endsWith('/')) {
        meta_base = "${meta_base}/"
    }

    Channel.from(params.input.metadata.files.collect { k, v ->
                           [k, file("${meta_base}${v}") ] })
                                             .into { METADATA_PRE; METADATA_PRE_VERBOSE }

    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA
                   .join(METADATA_PRE, remainder: true)
                   .into { METADATA_RIBO; METADATA_VERBOSE }
}
else {
    METADATA_PRE = Channel.from([null])
    METADATA_PRE_VERBOSE = Channel.from([null])

    GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_FINAL_FOR_METADATA
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

do_create_ribo = params.get('do_ribo_creation', false)

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

    Channel.from(params.get('rnaseq', [:]).fastq.collect { k, v ->
        v.collect { z -> [k, v.indexOf(z) + 1,
                                                   file("${rnaseq_fastq_base}${z}")] }  })
    .flatten().collate(3).into {  RNASEQ_FASTQ
                             RNASEQ_FASTQ_VERBOSE
                             RNASEQ_FASTQ_FASTQC
                             RNASEQ_FASTQ_CLIP
                             RNASEQ_FASTQ_EXISTENCE}

    if (params.do_check_file_existence) {

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
    set -o pipefail
    bowtie2 ${params.rnaseq.filter_arguments} \
            -x ${bowtie2_index_base} -q ${fastq} \
            --threads ${task.cpus} \
            --al-gz ${sample}.${index}.aligned.filter.fastq.gz \
            --un-gz ${sample}.${index}.unaligned.filter.fastq.gz \
                     2> ${sample}.${index}.filter.log \
            | samtools view -b - \
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

    // RNA-seq genome alignment
    do_rnaseq_genome_trim = params.rnaseq.containsKey('genome_read_trim') &&
                           params.rnaseq.genome_read_trim.get('enabled', true) &&
                           params.rnaseq.genome_read_trim.containsKey('length')

    if (do_rnaseq_genome_trim) {
        // Add trimming step before genome alignment
        process rnaseq_genome_read_trim {
            storeDir get_storedir('read_trimming', true) + '/' + params.output.individual_lane_directory

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

        // Split trimmed reads channel
        RNASEQ_GENOME_TRIMMED.into { RNASEQ_GENOME_TRIMMED_FOR_ALIGNMENT; RNASEQ_GENOME_TRIMMED_STORED }

        // Use trimmed reads for genome alignment
        RNASEQ_GENOME_INPUT_CHANNEL = RNASEQ_GENOME_TRIMMED_FOR_ALIGNMENT
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

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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
    set -o pipefail
    hisat2 ${params.get('rnaseq', [:]).get('hisat2_arguments', params.alignment_arguments.genome)} \\
           -x ${genome_base} -U ${fastq} \\
           -p ${task.cpus} \\
           --al-gz ${sample}.${index}.rnaseq_genome_alignment.aligned.fastq.gz \\
           --un-gz ${sample}.${index}.rnaseq_genome_alignment.unaligned.fastq.gz \\
                2> ${sample}.${index}.rnaseq_genome_alignment.log \\
           | samtools view -b - \\
           | samtools addreplacerg -r "ID:${sample}.${index}" -r "SM:${sample}" -r "PL:ILLUMINA" - \\
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

        """
    samtools view -bq ${params.mapping_quality_cutoff} ${bam} > ${sample}.${index}.rnaseq_genome_alignment.qpass.bam && \\
    samtools index ${sample}.${index}.rnaseq_genome_alignment.qpass.bam
    samtools view -c ${sample}.${index}.rnaseq_genome_alignment.qpass.bam > ${sample}.${index}.rnaseq_genome.qpass.count
    samtools flagstat ${sample}.${index}.rnaseq_genome_alignment.qpass.bam > ${sample}.${index}.rnaseq_genome_alignment.qpass.stats
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
    
    // Group RNA-seq genome alignment outputs for merging
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

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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
    (for f in ${aligned_fastq}; do if [ -s "\$f" ] && gzip -t "\$f" 2>/dev/null; then zcat "\$f"; fi; done) | gzip -c > ${sample}.rnaseq_genome.aligned.fastq.gz && \\
    (for f in ${unaligned_fastq}; do if [ -s "\$f" ] && gzip -t "\$f" 2>/dev/null; then zcat "\$f"; fi; done) | gzip -c > ${sample}.rnaseq_genome.unaligned.fastq.gz && \\
    cat ${alignment_log} > ${sample}.rnaseq_genome.log && \\
    rfc hisat2-log-to-csv \\
                       -l ${sample}.rnaseq_genome.log -n ${sample} -p rnaseq_genome \\
                       -o ${sample}.rnaseq_genome.csv
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
        samtools merge ${sample}.rnaseq_genome.qpass.merged.bam ${bam_files} && \
        samtools index ${sample}.rnaseq_genome.qpass.merged.bam
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
        storeDir get_storedir('alignment_ribo', true) + '/' + params.output.merged_lane_directory

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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
        storeDir  get_storedir('alignment_ribo', true) + '/' + params.output.individual_lane_directory

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

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

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
            -t ${genome_log} \
            -q ${qpass_count} \
            -d ${dedup_count} \
            --label-prefix genome \
            -o ${sample}.${index}.rnaseq_genome_individual.csv
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
        storeDir get_storedir('alignment_ribo', true) + '/' + params.output.merged_lane_directory

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
        python3 ${workflow.projectDir}/scripts/extract_reads_from_dedup_bed.py \\
            --bam ${qpass_bam} \\
            --bed ${post_dedup_bed} \\
            --output ${experiment}.rnaseq_genome.post_dedup.bam && \\
        samtools index ${experiment}.rnaseq_genome.post_dedup.bam
        """
    }

    process rnaseq_create_strand_specific_bigwigs {
        storeDir get_storedir('alignment_ribo', true) + '/bigwigs/' + params.output.merged_lane_directory

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_bigwig'

        input:
        set val(experiment), file(bam), file(bai) from RNASEQ_GENOME_DEDUP_BAM_FOR_BIGWIG
            .join(RNASEQ_GENOME_DEDUP_BAI_FOR_BIGWIG)

        output:
        set val(experiment), file("${experiment}.rnaseq.dedup.plus.bigWig"), \
                            file("${experiment}.rnaseq.dedup.minus.bigWig") \
            into RNASEQ_STRAND_SPECIFIC_BIGWIGS

        when:
        rnaseq_dedup_method == 'position'

        """
        samtools index ${bam}

        if [ "${rnaseq_library_strandedness}" == "forward" ]; then
            bamCoverage -b ${bam} -o ${experiment}.rnaseq.dedup.plus.bigWig \
                --filterRNAstrand reverse --binSize 1 -p ${task.cpus}
            bamCoverage -b ${bam} -o ${experiment}.rnaseq.dedup.minus.bigWig \
                --filterRNAstrand forward --binSize 1 -p ${task.cpus}
        else
            bamCoverage -b ${bam} -o ${experiment}.rnaseq.dedup.plus.bigWig \
                --filterRNAstrand forward --binSize 1 -p ${task.cpus}
            bamCoverage -b ${bam} -o ${experiment}.rnaseq.dedup.minus.bigWig \
                --filterRNAstrand reverse --binSize 1 -p ${task.cpus}
        fi
        """
    }

///////////////////////////////////////////////////////////////////////////////
//  E N D   R N A - S E Q   B I G W I G   G E N E R A T I O N  (MERGED)
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  R N A - S E Q   B I G W I G   G E N E R A T I O N  (NO DEDUP)
///////////////////////////////////////////////////////////////////////////////

    // Create bigWig files directly from quality-filtered BAMs when dedup_method == 'none'
    process rnaseq_create_nodedup_bigwigs {
        storeDir get_storedir('alignment_ribo', true) + '/bigwigs/' + params.output.merged_lane_directory

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_bigwig'

        input:
        set val(sample), file(bam), file(bai) from RNASEQ_GENOME_FOR_NODEDUP_BIGWIG

        output:
        set val(sample), file("${sample}.rnaseq.nodedup.plus.bigWig"), \
                         file("${sample}.rnaseq.nodedup.minus.bigWig") \
            into RNASEQ_NODEDUP_STRAND_SPECIFIC_BIGWIGS

        when:
        rnaseq_dedup_method == 'none'

        """
        samtools index ${bam}

        if [ "${rnaseq_library_strandedness}" == "forward" ]; then
            bamCoverage -b ${bam} -o ${sample}.rnaseq.nodedup.plus.bigWig \
                --filterRNAstrand reverse --binSize 1 -p ${task.cpus}
            bamCoverage -b ${bam} -o ${sample}.rnaseq.nodedup.minus.bigWig \
                --filterRNAstrand forward --binSize 1 -p ${task.cpus}
        else
            bamCoverage -b ${bam} -o ${sample}.rnaseq.nodedup.plus.bigWig \
                --filterRNAstrand forward --binSize 1 -p ${task.cpus}
            bamCoverage -b ${bam} -o ${sample}.rnaseq.nodedup.minus.bigWig \
                --filterRNAstrand reverse --binSize 1 -p ${task.cpus}
        fi
        """
    }

///////////////////////////////////////////////////////////////////////////////
//  E N D   R N A - S E Q   B I G W I G   G E N E R A T I O N  (NO DEDUP)
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
        publishDir get_rnaseq_publishdir("stats"), mode: 'copy'

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

    input:
        file(stat_table) from RNASEQ_GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

    output:
        file('rnaseq_individual_stats.csv') \
          into COMBINED_INDIVIDUAL_RNASEQ_GENOME_ALIGNMENT_STATS

        script:
        if (stat_table.size() == 0) {
            """
            echo "No RNA-seq individual statistics data available" > rnaseq_individual_stats.csv
            echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> rnaseq_individual_stats.csv
            """
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

    // Prepare input channel - normalize structure to always have 3 elements
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

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

        input:
        set val(sample), file(stat_files), val(merged_dedup_count_file) from RNASEQ_GENOME_STATS_INPUT

        output:
        set val(sample), file("${sample}.rnaseq_genome_merged.csv") into RNASEQ_GENOME_MERGED_ALIGNMENT_STATS

        script:
        // Check if merged_dedup_count_file exists and is not empty
        merged_count_file_exists = merged_dedup_count_file != null && merged_dedup_count_file.toString() != 'null' && merged_dedup_count_file.toString() != ''
        
        if (rnaseq_dedup_method == 'position' && merged_count_file_exists) {
            """
            # Sum individual stats first
            rfc sum-stats -n ${sample} -o ${sample}.rnaseq_genome_merged.tmp.csv ${stat_files}
            
            # Update the merged dedup count from the actual merged BAM count
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
        publishDir get_rnaseq_publishdir("stats"), mode: 'copy'

        beforeScript 'eval "$(conda shell.bash hook)" && conda activate ribo_genome'

    input:
        file(stat_table) from RNASEQ_GENOME_MERGED_ALIGNMENT_STATS_COLLECTED

    output:
        file('rnaseq_stats.csv') \
          into COMBINED_MERGED_RNASEQ_GENOME_ALIGNMENT_STATS

        script:
        if (stat_table.size() == 0) {
            // Create empty stats file when no input is available
            """
            echo "No RNA-seq merged statistics data available" > rnaseq_stats.csv
            echo "sample,total_reads,clipped_reads,filtered_out,filter_kept,genome_aligned_once,genome_aligned_many,genome_total_aligned,genome_unaligned,genome_qpass_aligned_reads,genome_after_dedup" >> rnaseq_stats.csv
            """
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
