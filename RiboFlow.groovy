/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

/*
Developed and tested on:
N E X T F L O W  ~  version 19.04.1
*/

// Default CPU setting - processes will use this unless overridden

// TODO
// DEDUPLICATE IN RNA-SEQ
// UPDATE DOCKER IMAGE
// Simplify Output

////////////////////////////////////////////////////////////////////////////////
////// General Function Definitions ////////////////////////////////////////////

String get_storedir(output_type){
    new File( params.output.intermediates.base,
              params.output.intermediates.get(output_type, output_type) )
							.getCanonicalPath()
}

String get_publishdir(output_type){
    new File( params.output.output.base,
              params.output.output.get(output_type, output_type) )
							.getCanonicalPath()
}

String get_dedup_method(String dedup_arg, String dedup_old){
  /*
  dedup_arg: User provided dedup method
  dedup_old: In the first version, dedup vas a boolean parameter
             We provide back-compatibility in the absenceo of dedup_arg

  Possible dedup_arg types are
    umi_tools: deduplication is based on umi umi_tools
               in this case looks for UMI_EXTRACT and UMI_DEDUP
               parameters if available to use them.

    position: position based deduplication.
              This is the method in the first version of ReiboFlow

    none:     Do not deduplicate

  */
  def valid_methods = ["position", "umi_tools", "none"]

  dedup_param = dedup_arg.toLowerCase()

  if(dedup_param != "none"){

    if(dedup_param in valid_methods){
      return(dedup_param)
    }
    else {
      println("Invalid deduplication method " + dedup_param + " . Valid methods are: ")
      println( valid_methods.join(",") )
      System.exit(1)
    }
  }
  else {
    // We make it compatible with the earlier versions.
    // So if the dedup flag is true, and dedup_method is not provided
    // it should be
    // position based dedup
    if( dedup_old.toLowerCase() != "false" ){
      return("position")
    }
    else{
      return("none")
    }
  }
}

////// General Function Definitions ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

dedup_method = get_dedup_method(params.get("dedup_method", "none").toString(),
                                params.get("deduplicate", false).toString() )

// to be deleted
//println(dedup_method)
//println(params.umi_tools_extract_arguments)
//System.exit(0)


fastq_base = params.input.get("fastq_base", "")
if(! fastq_base.endsWith("/") && fastq_base != ""){
	fastq_base = "${fastq_base}/"
}

// Group input files into a list of tuples where each item is
// [ sample, fileindex, path_to_fastq_file]
Channel.from(params.input.fastq.collect{k,v ->
	              v.collect{ z -> [k, v.indexOf(z) + 1,
									               file("${fastq_base}${z}")] }  })
	.flatten().collate(3).into{  INPUT_SAMPLES_VERBOSE;
		                           INPUT_SAMPLES_MD5;
	                             INPUT_SAMPLES_EXISTENCE;
	                             INPUT_SAMPLES_FASTQC;
															 INPUT_SAMPLES_CLIP;
	                             INPUT_SAMPLES_LOG;
															 INPUT_SAMPLES_READ_LENGTH;
														   INPUT_FOR_METADATA}



// Create a log file of index <-> fastq-file correspondence

INPUT_SAMPLES_LOG.flatMap{ sample, index, fastq -> "${sample}\t${index}\t${fastq}" }
    .collectFile(name: 'correspondence.txt', newLine: true)
    .set{INPUT_SAMPLES_LOG_FILES}

// Move the above correspondence file to an output folder via a process
process write_fastq_correspondence{

    executor 'local'

	publishDir get_publishdir("stats"), mode: 'move'

	input:
	file(correspondence) from INPUT_SAMPLES_LOG_FILES

	output:
	file("index_fastq_correspondence.txt")

	"""
    cat ${correspondence} > index_fastq_correspondence.txt
	"""
}


////////////////////////////////////////////////////////////////////////////////
////// Check File Existence ////////////////////////////////////////////////////

boolean file_exists(file_path) {
    this_file = file(file_path)
    assert this_file.exists()
    return true
}

boolean hisat2_ref_exists(hisat2_ref) {
    Channel.from( ["1.ht2", "2.ht2", "3.ht2","4.ht2","5.ht2","6.ht2", "7.ht2", "8.ht2"])
    .map{ this_suffix -> file_exists( "${hisat2_ref}.${this_suffix}".replaceAll('\\*', "") ) }
    return true
}

boolean bt2_ref_exists(bt2_ref) {
    Channel.from( ["1.bt2", "2.bt2", "3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    .map{ this_suffix -> file_exists( "${bt2_ref}.${this_suffix}".replaceAll('\\*', "") ) }
    return true
}

if(params.do_check_file_existence){
   // Make Sure Fastq Files Exist
   INPUT_SAMPLES_EXISTENCE.map{ sample, index, this_file -> file_exists(this_file) }

   // Make Sure bt2 and hisat reference files exist.
   bt2_ref_exists( params.input.reference.filter )
   if( params.input.reference.containsKey("transcriptome") ){
       bt2_ref_exists( params.input.reference.transcriptome )
   }
   if( params.input.reference.get("genome", false) ){
       hisat2_ref_exists( params.input.reference.genome )
   }

   if( params.input.reference.get("post_genome", false) ){
       bt2_ref_exists( params.input.reference.post_genome )
   }

   file_exists(params.input.reference.regions)
   file_exists(params.input.reference.transcript_lengths)

   root_meta_file = params.input.get("root_meta", false)
   if( root_meta_file ){
     file_exists(root_meta_file)
   }
}

////// Check File Existence ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////     P R O C E S S E S     /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* RAW_FASTQC */

process raw_fastqc{

	publishDir get_publishdir("fastqc")+"/raw", mode: 'copy'

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

// RAW_FASTQC
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* CLIP */

process clip{

    storeDir get_storedir("clip")

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

// CLIP
///////////////////////////////////////////////////////////////////////////////////////

CLIP_OUT.into{ CLIP_OUT_FASTQC; CLIP_OUT_DEDUP; CLIP_OUT_FILTER ; CLIP_OUT_READ_LENGTH}


///////////////////////////////////////////////////////////////////////////////////////
/* EXTRACT UMI via UMI_tools */

process extract_umi_via_umi_tools{
  storeDir get_storedir("umi_tools") + "/" + params.output.merged_lane_directory

  input:
  set val(sample), val(index), file(fastq) from CLIP_OUT_DEDUP

  output:
  set val(sample), val(index), file("${sample}.${index}.umi_extracted.fastq.gz") into UMI_EXTRACT_OUT
  set val(sample), val(index), file("${sample}.${index}.umi_extracted.log") into UMI_EXTRACT_LOG

  when:
  dedup_method == "umi_tools"

  """
  umi_tools extract -I ${fastq} -S ${sample}.${index}.umi_extracted.fastq.gz \
     -L ${sample}.${index}.umi_extracted.log \
     ${params.get("umi_tools_extract_arguments", "")}
  """
}

/* EXTRACT UMI via UMI_tools */
///////////////////////////////////////////////////////////////////////////////////////



if(dedup_method == "none" || dedup_method == "position"){
  CLIP_OUT_FILTER.set{ FILTER_INPUT_FASTQ }
}
else if(dedup_method == "umi_tools"){
  UMI_EXTRACT_OUT.set{ FILTER_INPUT_FASTQ }
}

///////////////////////////////////////////////////////////////////////////////////////
/* CLIPPED FASTQC */

process clipped_fastqc{

    publishDir get_publishdir("fastqc") + "/clipped", mode: 'copy'

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

// CLIPPED FASTQC
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* FILTER */

// Reads are mapped against (typically) rRNA, tRNA, adapter sequneces and etc
// that have no use for downstream analysis
// So we take the UNaligned reads from this process and use it for downstream processing

FILTER_INDEX = Channel.from([[
             params.input.reference.filter
                .split('/')[-1]
                .replaceAll('\\*$', "")
                .replaceAll('\\.$', ""),
             file(params.input.reference.filter),
            ]])

process filter{

	storeDir get_storedir("filter")

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

FILTER_ALIGNED.into{FILTER_ALIGNED_FASTQ_READ_LENGTH;
                    FILTER_ALIGNED_FASTQ_FASTQC}

// Check if transcriptome alignment should be performed
do_align_transcriptome = params.input.reference.containsKey("transcriptome")

if(do_align_transcriptome){
    FILTER_UNALIGNED.into{FILTER_UNALIGNED_FASTQ_READ_LENGTH;
                          FILTER_UNALIGNED_FASTQ_FASTQC;
                          FILTER_UNALIGNED_TRANSCRIPTOME;
                          FILTER_UNALIGNED_GENOME}
} else {
    FILTER_UNALIGNED.into{FILTER_UNALIGNED_FASTQ_READ_LENGTH;
                          FILTER_UNALIGNED_FASTQ_FASTQC;
                          FILTER_UNALIGNED_GENOME}
}

// FILTER
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
/* TRANSCRIPTOME ALIGNMENT */

if(do_align_transcriptome){

TRANSCRIPTOME_INDEX = Channel.from([[
             params.input.reference.transcriptome
                .split('/')[-1]
                .replaceAll('\\*$', "")
                .replaceAll('\\.$', ""),
             file(params.input.reference.transcriptome),
            ]])


process transcriptome_alignment{

    storeDir get_storedir("transcriptome_alignment") + "/" + params.output.individual_lane_directory

    input:
    set val(sample), val(index), file(fastq) from FILTER_UNALIGNED_TRANSCRIPTOME
		set val(transcriptome_reference), file(transcriptome_Reference_files) \
		            from TRANSCRIPTOME_INDEX.first()

    output:
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.bam") \
        into TRANSCRIPTOME_ALIGNMENT_BAM_PRE
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.bam.bai") \
        into TRANSCRIPTOME_ALIGNMENT_BAI
    set val(sample), val(index), file("${sample}.${index}.aligned.transcriptome_alignment.fastq.gz") \
        into TRANSCRIPTOME_ALIGNMENT_ALIGNED
    set val(sample), val(index), file("${sample}.${index}.unaligned.transcriptome_alignment.fastq.gz") \
        into TRANSCRIPTOME_ALIGNMENT_UNALIGNED
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.log") \
        into TRANSCRIPTOME_ALIGNMENT_LOG
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.stats") \
        into TRANSCRIPTOME_ALIGNMENT_STATS

    """
    bowtie2 ${params.alignment_arguments.transcriptome} \
            -x ${transcriptome_reference} -q ${fastq} \
            --threads ${task.cpus} \
            --al-gz ${sample}.${index}.aligned.transcriptome_alignment.fastq.gz \
            --un-gz ${sample}.${index}.unaligned.transcriptome_alignment.fastq.gz \
                     2> ${sample}.${index}.transcriptome_alignment.log \
            | samtools view -bS -o ${sample}.${index}.tmp.transcriptome_alignment.bam \
            && samtools addreplacerg -r ID:${sample}.${index} -@ ${task.cpus} \
                                     ${sample}.${index}.tmp.transcriptome_alignment.bam \
            | samtools sort -@ ${task.cpus} -o ${sample}.${index}.transcriptome_alignment.bam \
            && samtools index -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.bam \
            && samtools idxstats -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.bam  > \
               ${sample}.${index}.transcriptome_alignment.stats

    rm ${sample}.${index}.tmp.transcriptome_alignment.bam
    """
}

// TRANSCRIPTOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////////////

TRANSCRIPTOME_ALIGNMENT_BAM_PRE.into{ TRANSCRIPTOME_ALIGNMENT_BAM;
	                                    TRANSCRIPTOME_ALIGNMENT_BAM_MERGE;
                                      TRANSCRIPTOME_ALIGNMENT_BAM_FOR_QUALITY}

///////////////////////////////////////////////////////////////////////////////////////
/* QUALITY FILTER */

process quality_filter{

	storeDir get_storedir("quality_filter")

	input:
	set val(sample), val(index), file(bam) from TRANSCRIPTOME_ALIGNMENT_BAM_FOR_QUALITY

	output:
	set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.qpass.bam") \
        into TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_PRE
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.qpass.bam.bai") \
        into TRANSCRIPTOME_ALIGNMENT_QPASS_BAI
    set val(sample), val(index), file("${sample}.${index}.qpass.count") \
        into TRANSCRIPTOME_QPASS_COUNTS
    set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.qpass.stats") \
        into TRANSCRIPTOME_ALIGNMENT_QPASS_STATS

	"""
	samtools view -b -q ${params.mapping_quality_cutoff} ${bam}\
	| samtools sort -@ ${task.cpus} -o ${sample}.${index}.transcriptome_alignment.qpass.bam \
	&& samtools view -b -c ${sample}.${index}.transcriptome_alignment.qpass.bam > ${sample}.${index}.qpass.count \
	&& samtools index -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.qpass.bam \
	&& samtools idxstats -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.qpass.bam  > \
               ${sample}.${index}.transcriptome_alignment.qpass.stats
	"""
}

TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_PRE.into{ QPASS_BAM_READ_LENGTH;
	                                          TRANSCRIPTOME_ALIGNMENT_QPASS_BAM;
                                            TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_MERGE;
                                            TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_DEDUP;
                                            TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_ARRANGE}

// We need to copy output channels of transcriptome alignment
// for merging and variaous steps of downstream processing

TRANSCRIPTOME_ALIGNMENT_BAI.into{ TRANSCRIPTOME_ALIGNMENT_BAI_MERGE ;
                                  TRANSCRIPTOME_ALIGNMENT_BAI_REGION_COUNT}

TRANSCRIPTOME_ALIGNMENT_ALIGNED.into{ TRANSCRIPTOME_ALIGNMENT_ALIGNED_MERGE ;
                                  TRANSCRIPTOME_ALIGNMENT_ALIGNED_LENGTH ;
                                  TRANSCRIPTOME_ALIGNMENT_ALIGNED_FASTQC }

TRANSCRIPTOME_ALIGNMENT_UNALIGNED.into{ TRANSCRIPTOME_ALIGNMENT_UNALIGNED_MERGE ;
	                              TRANSCRIPTOME_ALIGNMENT_UNALIGNED_GENOME ;
                                  TRANSCRIPTOME_ALIGNMENT_UNALIGNED_LENGTH ;
                                  TRANSCRIPTOME_ALIGNMENT_UNALIGNED_FASTQC }

TRANSCRIPTOME_ALIGNMENT_LOG.into{ TRANSCRIPTOME_ALIGNMENT_LOG_MERGE ;
                                  TRANSCRIPTOME_ALIGNMENT_LOG_TABLE  }

TRANSCRIPTOME_ALIGNMENT_STATS.into{ TRANSCRIPTOME_ALIGNMENT_STATS_MERGE ;
                                    TRANSCRIPTOME_ALIGNMENT_STATS_TABLE  }

TRANSCRIPTOME_QPASS_COUNTS.into{TRANSCRIPTOME_QPASS_COUNTS_FOR_INDEX;
	                            TRANSCRIPTOME_QPASS_COUNTS_FOR_TABLE}

///////////////////////////////////////////////////////////////////////////////////////
/* MERGE TRANSCRIPTOME ALIGNMENT */

TRANSCRIPTOME_ALIGNMENT_BAM_MERGE.map{sample, index, bam -> [sample, bam]}.groupTuple()
    .set{ TRANSCRIPTOME_ALIGNMENT_GROUPED_BAM }

TRANSCRIPTOME_ALIGNMENT_ALIGNED_MERGE.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
    .set{ TRANSCRIPTOME_ALIGNMENT_GROUPED_ALIGNED }

TRANSCRIPTOME_ALIGNMENT_UNALIGNED_MERGE.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
    .set{ TRANSCRIPTOME_ALIGNMENT_GROUPED_UNALIGNED }

TRANSCRIPTOME_ALIGNMENT_LOG_MERGE.map{sample, index, log -> [sample, log]}.groupTuple()
    .set{ TRANSCRIPTOME_ALIGNMENT_GROUPED_LOG }


TRANSCRIPTOME_ALIGNMENT_GROUPED_BAM.join(TRANSCRIPTOME_ALIGNMENT_GROUPED_ALIGNED)
                                   .join(TRANSCRIPTOME_ALIGNMENT_GROUPED_UNALIGNED)
                                   .join(TRANSCRIPTOME_ALIGNMENT_GROUPED_LOG)
                                   .into{TRANSCRIPTOME_ALIGNMENT_GROUPED_JOINT;
                                         TRANSCRIPTOME_ALIGNMENT_GROUPED_JOINED_VERBOSE}


process merge_transcriptome_alignment{

	storeDir get_storedir("transcriptome_alignment") + "/" + params.output.merged_lane_directory

  input:
  set val(sample), file(bam_list), file(aligned_fastq_list), file(unaligned_fastq_list), file(log_list) \
             from TRANSCRIPTOME_ALIGNMENT_GROUPED_JOINT

  output:
  file("${sample}.transcriptome.bam") into TRANSCRIPTOME_ALIGNMENT_MERGED_BAM
  file("${sample}.transcriptome.bam.bai") into TRANSCRIPTOME_ALIGNMENT_MERGED_BAI
  file("${sample}.transcriptome.aligned.fastq.gz") into TRANSCRIPTOME_ALIGNMENT_MERGED_ALIGNED
  file("${sample}.transcriptome.unaligned.fastq.gz") into TRANSCRIPTOME_ALIGNMENT_MERGED_UNALIGNED
  file("${sample}.transcriptome.log") into TRANSCRIPTOME_ALIGNMENT_MERGED_LOG

  """
  samtools merge -@ ${task.cpus} ${sample}.transcriptome.bam ${bam_list} && \
  samtools index ${sample}.transcriptome.bam
  cat ${aligned_fastq_list} > ${sample}.transcriptome.aligned.fastq.gz
  cat ${unaligned_fastq_list} > ${sample}.transcriptome.unaligned.fastq.gz
    rfc merge bowtie2-logs --out ${sample}.transcriptome.log ${log_list}
	"""

}

// MERGE TRANSCRIPTOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* TRANSCRIPTOME INDIVIDUAL  FASTQC */

process transcriptome_aligned_individual_fastqc{
    publishDir get_publishdir("fastqc") + "/transcriptome_aligned", mode: 'copy'

    input:
    set val(sample), val(index), file(fastq)  from TRANSCRIPTOME_ALIGNMENT_ALIGNED_FASTQC

    output:
    set val(sample), file("${sample}.${index}.transcriptome.aligned_fastqc.html"),
                       file("${sample}.${index}.transcriptome.aligned_fastqc.zip") \
                        into TRANSCRIPTOME_ALIGNED_INDIVIDUAL_FASTQC_OUT

    when:
    params.do_fastqc

    """
    if [ ! -f ${sample}.${index}.transcriptome.aligned.fastq.gz ]; then
       ln -s ${fastq} ${sample}.${index}.transcriptome.aligned.fastq.gz
    fi
    fastqc ${sample}.${index}.transcriptome.aligned.fastq.gz --outdir=\$PWD -t ${task.cpus}
    """

}

process transcriptome_unaligned_individual_fastqc{
    publishDir get_publishdir("fastqc") + "/transcriptome_unaligned", mode: 'copy'

    input:
    set val(sample), val(index), file(fastq)  from TRANSCRIPTOME_ALIGNMENT_UNALIGNED_FASTQC


    output:
    set val(sample), file("${sample}.${index}.transcriptome.unaligned_fastqc.html"),
                       file("${sample}.${index}.transcriptome.unaligned_fastqc.zip") \
                        into TRANSCRIPTOME_UNALIGNED_INDIVIDUAL_FASTQC_OUT
    when:
    params.do_fastqc

    """
    if [ ! -f ${sample}.${index}.transcriptome.unaligned.fastq.gz ]; then
       ln -s ${fastq} ${sample}.${index}.transcriptome.unaligned.fastq.gz
    fi
    fastqc ${sample}.${index}.transcriptome.unaligned.fastq.gz --outdir=\$PWD -t ${task.cpus}
    """

}

} // end of if(do_align_transcriptome)
else {
    // Create empty channels when transcriptome alignment is disabled
    TRANSCRIPTOME_ALIGNMENT_QPASS_BAM = Channel.empty()
    TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_MERGE = Channel.empty()
    TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_DEDUP = Channel.empty()
    TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_ARRANGE = Channel.empty()
    QPASS_BAM_READ_LENGTH = Channel.empty()
    TRANSCRIPTOME_QPASS_COUNTS_FOR_INDEX = Channel.empty()
    TRANSCRIPTOME_QPASS_COUNTS_FOR_TABLE = Channel.empty()
    TRANSCRIPTOME_ALIGNMENT_LOG_TABLE = Channel.empty()
    TRANSCRIPTOME_ALIGNMENT_UNALIGNED_GENOME = Channel.empty()
    COMBINED_INDIVIDUAL_ALIGNMENT_STATS = Channel.empty()
    COMBINED_MERGED_ALIGNMENT_STATS = Channel.empty()
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* GENOME ALIGNMENT */

do_align_genome = params.input.reference.get("genome", false)

if(do_align_genome){

// Set up the input channel for genome alignment
// Genome alignment should always use filter unaligned reads regardless of transcriptome settings
GENOME_INPUT_CHANNEL = FILTER_UNALIGNED_GENOME

GENOME_INDEX = Channel.from([[
                params.input.reference.genome
                .split('/')[-1]
                .replaceAll('\\*$', "")
                .replaceAll('\\.$', ""),
             file(params.input.reference.genome),
            ]])


process genome_alignment{

	storeDir get_storedir("genome_alignment") + "/" + params.output.individual_lane_directory

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
              ${sample}.${index}.genome_alignment.stats
    """

}

GENOME_ALIGNMENT_ALIGNED.into{ GENOME_ALIGNMENT_ALIGNED_FASTQ_READ_LENGTH;
                               GENOME_ALIGNMENT_ALIGNED_MERGE;
                               GENOME_ALIGNMENT_ALIGNED_FASTQ_FASTQC }

GENOME_ALIGNMENT_UNALIGNED.into{ GENOME_ALIGNMENT_UNALIGNED_FASTQ_READ_LENGTH;
                                 GENOME_ALIGNMENT_UNALIGNED_MERGE;
                                 GENOME_ALIGNMENT_UNALIGNED_FASTQ_FASTQC;
                                 FOR_POST_GENOME }



// GENOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* MERGE GENOME ALIGNMENT */
GENOME_ALIGNMENT_LOG.into{ GENOME_ALIGNMENT_LOG_MERGE; GENOME_ALIGNMENT_LOG_TABLE; GENOME_ALIGNMENT_LOG_STATS }

GENOME_ALIGNMENT_LOG_TABLE
    .map{ sample, index, genome_log -> [ [sample, index], genome_log ] }
    .set{GENOME_ALIGNMENT_LOG_TABLE_INDEXED}


// Split genome alignment BAM for individual bed files, quality filtering, and merging
GENOME_ALIGNMENT_BAM.into{
    GENOME_ALIGNMENT_BAM_FOR_BED;
    GENOME_ALIGNMENT_BAM_FOR_QPASS;
    GENOME_ALIGNMENT_BAM_FOR_MERGE
}

// Apply quality filtering to genome BAM files
process genome_quality_filter{
    storeDir get_storedir("quality_filter")

    input:
    set val(sample), val(index), file(bam) from GENOME_ALIGNMENT_BAM_FOR_QPASS

    output:
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.bam") into GENOME_ALIGNMENT_QPASS_BAM
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.bam.bai") into GENOME_ALIGNMENT_QPASS_BAI
    set val(sample), val(index), file("${sample}.${index}.genome.qpass.count") into GENOME_QPASS_COUNTS
    set val(sample), val(index), file("${sample}.${index}.genome_alignment.qpass.stats") into GENOME_QPASS_STATS

    when:
    do_align_genome

    """
    samtools view -bq ${params.mapping_quality_cutoff} ${bam} > ${sample}.${index}.genome_alignment.qpass.bam && \
    samtools index ${sample}.${index}.genome_alignment.qpass.bam
    samtools view -c ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome.qpass.count
    samtools flagstat ${sample}.${index}.genome_alignment.qpass.bam > ${sample}.${index}.genome_alignment.qpass.stats
    """
}

// Split genome qpass BAM channel for various downstream processes
GENOME_ALIGNMENT_QPASS_BAM.into{ GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE;
                                  GENOME_ALIGNMENT_QPASS_BAM_FOR_STATS }

// Group genome qpass BAM files for merging (similar to transcriptome)
GENOME_ALIGNMENT_QPASS_BAM_FOR_MERGE.map{sample, index, bam -> [sample, bam]}.groupTuple()
   .set{ GENOME_ALIGNMENT_QPASS_GROUPED_BAM }

// Merge genome qpass BAM files (similar to transcriptome merge_bam_post_qpass)
process merge_genome_bam_post_qpass{
  storeDir get_storedir("genome_alignment") + "/" + params.output.merged_lane_directory

  input:
  set val(sample), file(bam_files) from GENOME_ALIGNMENT_QPASS_GROUPED_BAM

  output:
  set val(sample), file("${sample}.genome.qpass.merged.bam"),\
                   file("${sample}.genome.qpass.merged.bam.bai") \
      into GENOME_MERGED_BAM_QPASS_PRE_DEDUP

  when:
	dedup_method == "umi_tools"

  """
  samtools merge ${sample}.genome.qpass.merged.bam ${bam_files} && samtools index ${sample}.genome.qpass.merged.bam
  """
}

// Convert individual genome BAM files to BED format
process individual_genome_bam_to_bed{
    storeDir get_storedir("bam_to_bed")

    input:
    set val(sample), val(index), file(bam) from GENOME_ALIGNMENT_BAM_FOR_BED

    output:
    set val(sample), val(index), file("${sample}.${index}.genome.bed") into GENOME_INDIVIDUAL_BED_PRE_SPLIT
    set val(sample), val(index), file("${sample}.${index}.genome_nodedup_count.txt") \
       into GENOME_INDIVIDUAL_NODEDUP_COUNT

    when:
    do_align_genome

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

// Split genome individual BED channel for multiple uses
GENOME_INDIVIDUAL_BED_PRE_SPLIT.into{
    GENOME_INDIVIDUAL_BED_FOR_INDEX;
    GENOME_INDIVIDUAL_BED_FOR_SPLITTING
}

// Add sample index column to individual genome bed files
process add_sample_index_col_to_genome_bed{
    storeDir get_storedir("bam_to_bed")

    input:
    set val(sample), val(index), file(bed) from GENOME_INDIVIDUAL_BED_FOR_INDEX

    output:
    set val(sample), file("${sample}.${index}.genome.with_sample_index.bed") \
         into GENOME_BED_FOR_DEDUP_INDEX_COL_ADDED

    when:
    do_align_genome

    """
    awk -v this_sample=${sample}.${index} \
     '{ print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"this_sample) }' ${bed} \
          > ${sample}.${index}.genome.with_sample_index.bed
    """
}

// Group for merging
GENOME_ALIGNMENT_BAM_FOR_MERGE.map{sample, index, bam -> [sample, bam]}.groupTuple()
    .set{ GENOME_ALIGNMENT_GROUPED_BAM }

GENOME_ALIGNMENT_ALIGNED_MERGE.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
    .set{ GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ }

GENOME_ALIGNMENT_UNALIGNED_MERGE.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
    .set{ GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ }
GENOME_ALIGNMENT_LOG_MERGE.map{sample, index, log -> [sample, log]}.groupTuple()
    .set{ GENOME_ALIGNMENT_GROUPED_LOG }


GENOME_ALIGNMENT_GROUPED_BAM.join( GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ )
                            .join(GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ)
                            .join(GENOME_ALIGNMENT_GROUPED_LOG)
                            .set{ GENOME_ALIGNMENT_GROUPED_JOINT }


process merge_genome_alignment{

	storeDir get_storedir("genome_alignment") + "/" + params.output.merged_lane_directory

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

	"""
	samtools merge ${sample}.genome.bam ${bam} && samtools index ${sample}.genome.bam && \
    zcat ${aligned_fastq} | gzip -c > ${sample}.genome.aligned.fastq.gz && \
    zcat ${unaligned_fastq} | gzip -c > ${sample}.genome.unaligned.fastq.gz && \
    python3 ${workflow.projectDir}/merge-hisat2-logs.py -o ${sample}.genome.log ${alignment_log}
	"""

}

///////////////////////////////////////////////////////////////////////////////
//          G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////

// Set genome BAM channel for bed conversion
GENOME_ALIGNMENT_MERGED_BAM.set{ GENOME_BAM_FOR_DEDUP }

process genome_bam_to_bed{

	storeDir get_storedir("bam_to_bed")

	input:
	set val(sample), file(bam) from GENOME_BAM_FOR_DEDUP

	output:
	set val(sample), file("${sample}.genome.bed") into GENOME_BED_FOR_DEDUP_PRE
	set val(sample), file("${sample}.genome_nodedup_count.txt") \
	   into GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP

	when:
	dedup_method != "none"

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

GENOME_BED_FOR_DEDUP_PRE.into{
    GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION;
    GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_UMI
}

// Position-based genome deduplication
process genome_deduplicate_position{

	storeDir get_storedir("alignment_ribo")

	input:
	set val(sample), file(bed) from GENOME_BED_FOR_DEDUP_MERGED_PRE_DEDUP_POSITION

	output:
	set val(sample), file("${sample}.genome.post_dedup.bed") \
	     into GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION

	when:
	dedup_method == "position"

	"""
	rfc dedup -i ${bed} -o ${sample}.genome.post_dedup.bed
	"""
}

// UMI-based genome deduplication requires BAM input

process genome_deduplicate_umi_tools{
  storeDir get_storedir("umi_tools") + "/" + params.output.merged_lane_directory

  input:
	set val(sample), file(bam), file(bai) from GENOME_MERGED_BAM_QPASS_PRE_DEDUP

  output:
  set val(sample), file("${sample}.genome.dedup.bam") into GENOME_UMI_TOOLS_DEDUP_BAM
  set val(sample), file("${sample}.genome.dedup.bed") into GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI
  set val(sample), file("${sample}.genome.dedup.log") into GENOME_UMI_TOOLS_DEDUP_LOG
  set val(sample), file("${sample}.genome.dedup.stats_edit_distance.tsv"),
                   file("${sample}.genome.dedup.stats_per_umi_per_position.tsv"),
                   file("${sample}.genome.dedup.stats_per_umi.tsv") \
                        into GENOME_UMI_TOOLS_DEDUP_STATS

  when:
  dedup_method == "umi_tools"

  """
  umi_tools dedup ${params.get("umi_tools_dedup_arguments", "")} \
              -I ${bam} --output-stats=${sample}.genome.dedup.stats -S ${sample}.genome.dedup.bam -L ${sample}.genome.dedup.log

  bamToBed -i ${sample}.genome.dedup.bam > ${sample}.genome.dedup.bed
  """
}

// Combine position and UMI deduplication outputs
GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_POSITION.mix(GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP_UMI)
    .set{GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP}

// Create individual genome post_dedup files by splitting merged BAM
GENOME_UMI_TOOLS_DEDUP_BAM
    .cross(GENOME_INDIVIDUAL_BED_FOR_SPLITTING.map{ sample, index, bed -> [sample, index] }.unique())
    .map{ merged_data, individual_data ->
        [individual_data[0], individual_data[1], merged_data[1]]
    }
    .set{ GENOME_DEDUP_BAM_FOR_SPLITTING }

process split_genome_dedup_bam_to_individual{
    storeDir get_storedir("alignment_ribo")

    input:
    set val(sample), val(index), file(merged_bam) from GENOME_DEDUP_BAM_FOR_SPLITTING

    output:
    set val(sample), val(index), file("${sample}.${index}.genome.post_dedup.bam") into GENOME_INDIVIDUAL_POST_DEDUP_BAM
    set val(sample), val(index), file("${sample}.${index}.genome.count_after_dedup.txt") into GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING

    when:
    dedup_method == "umi_tools" && do_align_genome

    """
    samtools view -B -r ${sample}.${index} ${merged_bam} -o ${sample}.${index}.genome.post_dedup.bam
    samtools view -c ${sample}.${index}.genome.post_dedup.bam > ${sample}.${index}.genome.count_after_dedup.txt
    """
}

// Create count file for deduplicated genome reads
process genome_count_deduplicated_reads{
    storeDir get_storedir("alignment_ribo")

    input:
    set val(sample), file(bed) from GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP

    output:
    set val(sample), file("${sample}.genome.dedup_count.txt") \
        into GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP

    when:
    dedup_method != "none"

    """
    wc -l ${bed} > ${sample}.genome.dedup_count.txt
    """
}

// Select the appropriate count channel based on deduplication method
if(dedup_method == "umi_tools"){
  // Use individual dedup counts from BAM splitting (has sample, index, count structure)
  GENOME_INDIVIDUAL_DEDUP_COUNT_FROM_SPLITTING.set{GENOME_INDIVIDUAL_DEDUP_COUNT}
}
else if(dedup_method != "none"){
  // Use merged dedup counts and convert to individual structure with dummy index
  GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP
    .map{ sample, count -> [sample, "1", count] }
    .set{GENOME_INDIVIDUAL_DEDUP_COUNT}
}
else{
  // Use merged no-dedup counts and convert to individual structure with dummy index
  GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP
    .map{ sample, count -> [sample, "1", count] }
    .set{GENOME_INDIVIDUAL_DEDUP_COUNT}
}

///////////////////////////////////////////////////////////////////////////////
//          E N D   G E N O M E   D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////



// MERGE  GENOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////

} // end of if(do_align_genome){

// END OF GENOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* BAM TO BED */


/*
We assume that duplicates are coming from PCR. So for each sample,
we merge the sequencing lanes ,
Add a sample.index column to the bed file.
sort the entire bed file by the first 3 columns ( sort -k1,1 -k2,2n -k3,3n )
Then we can deduplicate the entire file
Then separate the file based on the additional column that we added
*/

process bam_to_bed{

	storeDir get_storedir("bam_to_bed")

	when:
	do_align_transcriptome

	input:
	set val(sample), val(index), file(bam) from TRANSCRIPTOME_ALIGNMENT_QPASS_BAM

	output:
	set val(sample), val(index), file("${sample}.${index}.transcriptome.bed") into BAM_TO_BED
	set val(sample), val(index), file("${sample}.${index}.transcriptome_nodedup_count.txt") \
	   into INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP

    """
    if [ `samtools view -c ${bam}` -eq 0 ];
    then
       touch ${sample}.${index}.transcriptome.bed
    else
        bamToBed -i ${bam} > ${sample}.${index}.transcriptome.bed
    fi

    wc -l ${sample}.${index}.transcriptome.bed > ${sample}.${index}.transcriptome_nodedup_count.txt
    """

}

BAM_TO_BED.into{ BED_NODEDUP; BED_FOR_DEDUP; BED_FOR_INDEX_SEP_PRE }


process add_sample_index_col_to_bed{

	storeDir get_storedir("bam_to_bed")

	input:
    set val(sample), val(index), file(bed) from BED_FOR_DEDUP

	output:
	set val(sample), file("${sample}.${index}.with_sample_index.bed")\
	     into BED_FOR_DEDUP_INDEX_COL_ADDED

	"""
	awk -v newcol=${sample}.${index} '{print(\$0"\\t"newcol)}' ${bed}\
	   > ${sample}.${index}.with_sample_index.bed
	"""
}

BED_FOR_DEDUP_INDEX_COL_ADDED.groupTuple()
    .set{ BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED }

process merge_bed{

	storeDir get_storedir("bam_to_bed")

	input:
	set val(sample), file(bed_files) from BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED

	output:
	set val(sample), file("${sample}.merged.pre_dedup.bed") \
	    into MERGE_BED_OUT


  when:
	dedup_method == "position"

	"""
	cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.merged.pre_dedup.bed
	"""
}

MERGE_BED_OUT.into{BED_FOR_DEDUP_MERGED_PRE_DEDUP;
                   MERGED_BED_FOR_RIBO}


TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_MERGE.map{sample, index, bam -> [sample, bam]}.groupTuple()
   .set{ TRANSCRIPTOME_ALIGNMENT_QPASS_GROUPED_BAM }


process merge_bam_post_qpass{
  storeDir get_storedir("transcriptome_alignment") + "/" + params.output.merged_lane_directory

  input:
  set val(sample), file(bam_files) from TRANSCRIPTOME_ALIGNMENT_QPASS_GROUPED_BAM

  output:
  set val(sample), file("${sample}.merged.pre_dedup.bam"),\
                   file("${sample}.merged.pre_dedup.bam.bai") \
      into MERGED_BAM_PRE_DEDUP

  when:
	dedup_method == "umi_tools"

  """
  samtools merge ${sample}.merged.pre_dedup.bam ${bam_files}

  samtools index ${sample}.merged.pre_dedup.bam ${sample}.merged.pre_dedup.bam.bai
  """
}

///////////////////////////////////////////////////////////////////////////////
//          D E D U P L I C A T I O N

process deduplicate_position{

	storeDir  get_storedir("alignment_ribo")

	input:
	set val(sample), file(bed) from BED_FOR_DEDUP_MERGED_PRE_DEDUP

	output:
	set val(sample), file("${sample}.merged.post_dedup.bed") \
	     into BED_FOR_DEDUP_MERGED_POST_DEDUP

	when:
	dedup_method == "position"

	"""
	rfc dedup -i ${bed} -o ${sample}.merged.post_dedup.bed
	"""
}


process deduplicate_umi_tools{
  storeDir  get_storedir("umi_tools") + "/" + params.output.merged_lane_directory

  input:
	set val(sample), file(bam), file(bai) from MERGED_BAM_PRE_DEDUP

  output:
  set val(sample), file("${sample}.transcriptome.dedup.bam") into UMI_TOOLS_DEDUP_BAM

  set val(sample), file("${sample}.transcriptome.dedup.bed") into  UMI_TOOLS_BED_FOR_DEDUP_MERGED_POST_DEDUP

  set val(sample), file("${sample}.transcriptome.dedup.log") into UMI_TOOLS_DEDUP_LOG

  set val(sample), file("${sample}.transcriptome.dedup.stats_edit_distance.tsv"),
                   file("${sample}.transcriptome.dedup.stats_per_umi_per_position.tsv"),
                   file("${sample}.transcriptome.dedup.stats_per_umi.tsv") \
                        into UMI_TOOLS_DEDUP_STATS

  when:
  dedup_method == "umi_tools"

  """
  umi_tools dedup ${params.get("umi_tools_dedup_arguments", "")} \
              -I ${bam} --output-stats=${sample}.transcriptome.dedup.stats -S ${sample}.transcriptome.dedup.bam -L ${sample}.transcriptome.dedup.log

  bamToBed -i ${sample}.transcriptome.dedup.bam > ${sample}.transcriptome.dedup.bed
  """
}

//          D E D U P L I C A T I O N
///////////////////////////////////////////////////////////////////////////////

BED_FOR_DEDUP_MERGED_POST_DEDUP.into{BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_SEP;
                                     BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_RIBO}

BED_FOR_INDEX_SEP_PRE.map{ sample,index,file -> [sample, index] }
    .combine(BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_SEP, by:0)
    .set{ BED_FOR_INDEX_SEP_POST_DEDUP }


process separate_bed_post_dedup{

	storeDir  get_storedir("alignment_ribo")

	input:
	set val(sample), val(index), file(bed) from BED_FOR_INDEX_SEP_POST_DEDUP

	output:
	set val(sample), val(index), file("${sample}.${index}.transcriptome.post_dedup.bed") \
	   into BED_DEDUPLICATED
	set val(sample), val(index), file("${sample}.${index}.transcriptome.count_after_dedup.txt")\
	   into INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP

  when:
  dedup_method == "position"

	"""
	awk -v this_sample=${sample}.${index} \
	 '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} \
        > ${sample}.${index}.transcriptome.post_dedup.bed \
	  && wc -l ${sample}.${index}.transcriptome.post_dedup.bed > ${sample}.${index}.transcriptome.count_after_dedup.txt
	"""
}

TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_FOR_ARRANGE.map{ sample,index,file -> [sample, index] }
    .combine(UMI_TOOLS_DEDUP_BAM , by:0)
    .set{ BAM_FOR_INDEX_SEP_POST_DEDUP }

process separate_bam_post_dedup{
  storeDir  get_storedir("alignment_ribo")

  input:
  set val(sample), val(index), file(bam) from BAM_FOR_INDEX_SEP_POST_DEDUP

  output:
  set val(sample), val(index), file("${sample}.${index}.transcriptome.post_dedup.bam") \
	   into BAM_INDIVIDUAL_DEDUPLICATED
  set val(sample), val(index), file("${sample}.${index}.transcriptome.count_after_dedup.txt")\
  	 into BAM_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP


  when:
  dedup_method == "umi_tools"

  """
  samtools view -B -r ${sample}.${index} ${bam} -o ${sample}.${index}.transcriptome.post_dedup.bam

  samtools view -c ${sample}.${index}.transcriptome.post_dedup.bam > ${sample}.${index}.transcriptome.count_after_dedup.txt
  """
}


if(dedup_method == "position"){
  BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_RIBO
  .into{BED_FOR_SEPARATION; BED_FOR_RIBO; BED_FOR_RIBO_VERBOSE; BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT}

  // Create count file for deduplicated transcriptome reads
  process transcriptome_count_deduplicated_reads_position{
      storeDir get_storedir("alignment_ribo")

      input:
      set val(sample), file(bed) from BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT

      output:
      set val(sample), file("${sample}.transcriptome.dedup_count.txt") \
          into TRANSCRIPTOME_MERGED_DEDUP_COUNT

      """
      wc -l ${bed} > ${sample}.transcriptome.dedup_count.txt
      """
  }

  INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP
  .set{INDIVIDUAL_DEDUP_COUNT}
}
else if(dedup_method == "umi_tools"){
  UMI_TOOLS_BED_FOR_DEDUP_MERGED_POST_DEDUP
  .into{BED_FOR_SEPARATION; BED_FOR_RIBO; BED_FOR_RIBO_VERBOSE; BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT}

  // Create count file for deduplicated transcriptome reads
  process transcriptome_count_deduplicated_reads_umi{
      storeDir get_storedir("alignment_ribo")

      input:
      set val(sample), file(bed) from BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_COUNT

      output:
      set val(sample), file("${sample}.transcriptome.dedup_count.txt") \
          into TRANSCRIPTOME_MERGED_DEDUP_COUNT

      """
      wc -l ${bed} > ${sample}.transcriptome.dedup_count.txt
      """
  }

  BAM_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP
  .set{INDIVIDUAL_DEDUP_COUNT}
}
else if(dedup_method == "none"){
  MERGED_BED_FOR_RIBO
  .into{BED_FOR_SEPARATION; BED_FOR_RIBO; BED_FOR_RIBO_VERBOSE}

  INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP
  .set{INDIVIDUAL_DEDUP_COUNT}
}



///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
/* INDIVIDUAL ALIGNMENT STATS TABLE */

// Create indexed channels for clip and filter logs (needed for both transcriptome and genome stats)
CLIP_LOG.map{ sample, index, clip_log -> [ [sample, index], clip_log ] }
        .set{CLIP_LOG_INDEXED_TEMP}
FILTER_LOG.map{ sample, index, filter_log -> [ [sample, index], filter_log ] }
          .set{FILTER_LOG_INDEXED_TEMP}

// Always split channels for both uses - conditionals will determine if processes run
CLIP_LOG_INDEXED_TEMP.into{
    CLIP_LOG_INDEXED;
    CLIP_LOG_INDEXED_FOR_GENOME
}
FILTER_LOG_INDEXED_TEMP.into{
    FILTER_LOG_INDEXED;
    FILTER_LOG_INDEXED_FOR_GENOME
}

if(do_align_transcriptome){

// We need to group the log files by sample name and index
// than flatten that list and group again so that each
// entry can be emmited in groups of 6 for each task

TRANSCRIPTOME_ALIGNMENT_LOG_TABLE
    .map{ sample, index, transcriptome_log -> [ [sample, index], transcriptome_log ] }
    .set{TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED_TEMP}

TRANSCRIPTOME_QPASS_COUNTS_FOR_INDEX
    .map{ sample, index, qpass_count -> [ [sample, index], qpass_count ] }
    .set{TRANSCRIPTOME_QPASS_COUNTS_INDEXED}
INDIVIDUAL_DEDUP_COUNT
     .map{ sample, index, dedup_count -> [ [sample, index], dedup_count ] }
     .set{ INDIVIDUAL_DEDUP_COUNT_INDEXED }

// Only create transcriptome alignment log channels if transcriptome is enabled
    TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED_TEMP.into{
        TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED;
        TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED_FOR_GENOME
    }

} else {
    // Create empty channels when transcriptome is disabled
    TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED = Channel.empty()
    TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED_FOR_GENOME = Channel.empty()
    TRANSCRIPTOME_QPASS_COUNTS_INDEXED = Channel.empty()
    INDIVIDUAL_DEDUP_COUNT_INDEXED = Channel.empty()
}

// TRANSCRIPTOME STATS PROCESSING
if(do_align_transcriptome) {

CLIP_LOG_INDEXED.join(FILTER_LOG_INDEXED)
                .join(TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED)
                .join(TRANSCRIPTOME_QPASS_COUNTS_INDEXED)
                .join(INDIVIDUAL_DEDUP_COUNT_INDEXED)
                .flatten()
                .collate(7)
                .set{ INDIVIDUAL_ALIGNMENT_STATS_INPUT }

process individual_alignment_stats{
	/*
	Compiles statistics coming from the individual steps:
	cutadapt, filter, transcriptome and genome alignment,
	quality filtering and deduplication
	*/

	 executor 'local'

   storeDir get_storedir("stats") + "/transcriptome/individual"

   input:
   set val(sample), val(index), file(clip_log), file(filter_log),\
       file(transcriptome_log), file(qpass_count),\
       file(dedup_count)\
       from INDIVIDUAL_ALIGNMENT_STATS_INPUT

   output:
   set val(sample), val(index), file("${sample}.${index}.transcriptome_individual.csv") \
      into INDIVIDUAL_ALIGNMENT_STATS

   """
   rfc compile-step-stats \
	   -n ${sample}.${index} -c ${clip_log} \
     -f ${filter_log} -t ${transcriptome_log} \
		 -q ${qpass_count} \
     -d ${dedup_count} \
     -o ${sample}.${index}.transcriptome_individual.csv
   """

}

// INDIVIDUAL ALIGNMENT STATS
///////////////////////////////////////////////////////////////////////////////////////

INDIVIDUAL_ALIGNMENT_STATS
    .into{ INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION;
           INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING}

///////////////////////////////////////////////////////////////////////////////////////
/* COMBINE INDIVIDUAL ALIGNMENT STATS */

INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
    .map{ sample, index, stats_file -> stats_file }
    .toSortedList().set{INDIVIDUAL_ALIGNMENT_STATS_COLLECTED}

process combine_individual_alignment_stats{

  executor 'local'

	storeDir get_storedir("stats") + "/transcriptome/individual"

	input:
	file(stat_table) from INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

	output:
	file("transcriptome_individual_essential.csv") \
	      into COMBINED_INDIVIDUAL_ALIGNMENT_STATS

	"""
	  rfc merge overall-stats \
	   -o raw_combined_individual_aln_stats.csv \
	      ${stat_table} && \
    rfc stats-percentage \
	  -i raw_combined_individual_aln_stats.csv \
	  -o transcriptome_individual_essential.csv
	"""
}


// COMBINE INDIVIDUAL ALIGNMENT STATS
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* SUM INDIVIDUAL ALIGNMENT STATS */

/*
For each sample, sums up the stats coming from individual lanes
*/

INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING
    .map{ sample, index, file -> [ sample, file ] }
    .groupTuple()
    .into{ INDIVIDUAL_ALIGNMENT_STATS_GROUPED ;
           INDIVIDUAL_ALIGNMENT_STATS_GROUPED_VERBOSE }

process sum_individual_alignment_stats{

  executor 'local'

	storeDir get_storedir( "stats/transcriptome/merged" )

	input:
	set val(sample), file(stat_files) from INDIVIDUAL_ALIGNMENT_STATS_GROUPED

	output:
	set val(sample), file("${sample}.merged.alignment_stats.csv")\
	   into MERGED_ALIGNMENT_STATS

	"""
	rfc sum-stats -n ${sample}\
	  -o ${sample}.merged.alignment_stats.csv ${stat_files}
	"""
}

// SUM INDIVIDUAL ALIGNMENT STATS
////////////////////////////////////////////////////////////////////////////////

MERGED_ALIGNMENT_STATS.map{ sample, stats_file -> stats_file }
                      .toSortedList()
                      .set{ MERGED_ALIGNMENT_STATS_COLLECTED }

////////////////////////////////////////////////////////////////////////////////
/* COMBINE MERGED ALIGNMENT STATS */

process combine_merged_alignment_stats{

	storeDir get_storedir("stats") + "/transcriptome/merged"

	executor 'local'

	input:
	file(stat_files) from MERGED_ALIGNMENT_STATS_COLLECTED

	output:
	file("transcriptome_merged_essential.csv") into COMBINED_MERGED_ALIGNMENT_STATS

	"""
	rfc merge overall-stats \
	    -o raw_combined_merged_aln_stats.csv \
	    ${stat_files} && \
	rfc stats-percentage \
	  -i raw_combined_merged_aln_stats.csv \
	  -o transcriptome_merged_essential.csv
	"""
}

// COMBINE MERGED ALIGNMENT STATS
////////////////////////////////////////////////////////////////////////////////

} // end of if(do_align_transcriptome) block for stats processing

///////////////////////////////////////////////////////////////////////////////
/* GENOME STATS */

if(do_align_genome){

  // Create indexed channels for genome stats - following transcriptome pattern exactly
  GENOME_ALIGNMENT_LOG_STATS
      .map{ sample, index, log_file -> [ [sample, index], log_file] }
      .set{GENOME_ALIGNMENT_LOG_STATS_INDEXED}

  GENOME_QPASS_COUNTS
      .map{ sample, index, count_file -> [ [sample, index], count_file] }
      .set{GENOME_QPASS_COUNTS_INDEXED}

  GENOME_INDIVIDUAL_DEDUP_COUNT
      .map{ sample, index, count_file -> [ [sample, index], count_file] }
      .set{GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED}

  // Create genome alignment stats input channel
  CLIP_LOG_INDEXED_FOR_GENOME.join(FILTER_LOG_INDEXED_FOR_GENOME)
                .join(GENOME_ALIGNMENT_LOG_STATS_INDEXED)
                .join(GENOME_QPASS_COUNTS_INDEXED)
                .join(GENOME_INDIVIDUAL_DEDUP_COUNT_INDEXED)
                .flatten()
                .collate(7)
                .set{ GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT }

  // INDIVIDUAL GENOME ALIGNMENT STATS - Create equivalent of transcriptome stats
  process individual_genome_alignment_stats{
    storeDir get_storedir("stats") + "/genome/individual"

    input:
    set val(sample), val(index), file(clip_log), file(filter_log),\
        file(genome_log), file(qpass_count),\
        file(dedup_count)\
        from GENOME_INDIVIDUAL_ALIGNMENT_STATS_INPUT

    output:
    set val(sample), val(index), file("${sample}.${index}.genome_individual.csv") \
       into GENOME_INDIVIDUAL_ALIGNMENT_STATS

    """
    rfc compile-step-stats \
      -n ${sample}.${index} \
      -c ${clip_log} \
      -f ${filter_log} \
      -t ${genome_log} \
      -q ${qpass_count} \
      -d ${dedup_count} \
      -l genome \
      -o ${sample}.${index}.genome_individual.csv
    """
  }

  // Split genome individual stats
  GENOME_INDIVIDUAL_ALIGNMENT_STATS
      .into{ GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION;
             GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING}

  // Collect genome individual stats
  GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
      .map{ sample, index, stats_file -> stats_file }
      .toSortedList().set{GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED}

  // COLLECT INDIVIDUAL GENOME ALIGNMENT STATS
  process combine_individual_genome_alignment_stats{
    executor 'local'
    storeDir get_storedir("stats") + "/genome/individual"

    input:
    file(stat_table) from GENOME_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

    output:
    file("genome_individual_essential.csv") \
          into COMBINED_INDIVIDUAL_GENOME_ALIGNMENT_STATS

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

  // MERGED GENOME ALIGNMENT STATS - Follow transcriptome pattern of grouping individual stats
  GENOME_INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING
      .map{ sample, index, file -> [ sample, file ] }
      .groupTuple()
      .set{ GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED }

  process sum_individual_genome_alignment_stats{
    executor 'local'
    storeDir get_storedir("stats") + "/genome/merged"

    input:
    set val(sample), file(stat_files) from GENOME_INDIVIDUAL_ALIGNMENT_STATS_GROUPED

    output:
    set val(sample), file("${sample}.genome_merged.csv") into GENOME_MERGED_ALIGNMENT_STATS

    """
    rfc sum-stats -n ${sample}\
      -o ${sample}.genome_merged.csv ${stat_files}
    """
  }

  // Split and collect merged genome stats
  GENOME_MERGED_ALIGNMENT_STATS
      .into{ GENOME_MERGED_ALIGNMENT_STATS_FOR_COLLECTION;
             GENOME_MERGED_ALIGNMENT_STATS_FOR_GROUPING}

  GENOME_MERGED_ALIGNMENT_STATS_FOR_COLLECTION
      .map{ sample, stats_file -> stats_file }
      .toSortedList().set{GENOME_MERGED_ALIGNMENT_STATS_COLLECTED}

  // COMBINE MERGED GENOME ALIGNMENT STATS
  process combine_merged_genome_alignment_stats{
    executor 'local'
    storeDir get_storedir("stats") + "/genome/merged"

    input:
    file(stat_files) from GENOME_MERGED_ALIGNMENT_STATS_COLLECTED

    output:
    file("genome_merged_essential.csv") into COMBINED_MERGED_GENOME_ALIGNMENT_STATS

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

  // CREATE COMPREHENSIVE ALIGNMENT STATS (TRANSCRIPTOME + GENOME)
  // Combine transcriptome and genome essential merged stats into one comprehensive file
  process create_comprehensive_alignment_stats{
    executor 'local'
    storeDir get_storedir("log") + "/" + params.output.merged_lane_directory

    input:
    file(transcriptome_essential) from COMBINED_MERGED_ALIGNMENT_STATS
    file(genome_essential) from COMBINED_MERGED_GENOME_ALIGNMENT_STATS

    output:
    file("*.merged.alignment_stats.csv") into COMPREHENSIVE_MERGED_STATS

    script:
    """
    # Extract sample name from transcriptome file
    sample_name=\$(head -1 ${transcriptome_essential} | cut -d',' -f2)

    # Create comprehensive stats by combining transcriptome and genome data
    python3 << 'EOF'
import pandas as pd

# Read transcriptome stats
trans_df = pd.read_csv("${transcriptome_essential}", index_col=0)
sample_name = trans_df.columns[0]

# Read genome stats
genome_df = pd.read_csv("${genome_essential}", index_col=0)

# Combine both dataframes
combined_df = pd.concat([trans_df, genome_df])

# Remove duplicate rows (like total_reads, clipped_reads, etc that appear in both)
combined_df = combined_df[~combined_df.index.duplicated(keep='first')]

# Save comprehensive stats
combined_df.to_csv(f"{sample_name}.merged.alignment_stats.csv")
EOF
    """
  }

} // end of if(do_align_genome) for stats

///////////////////////////////////////////////////////////////////////////////
/* METADATA CHANNELS */
do_metadata = params.get("do_metadata", false) && params.input.get("metadata", false)

if( do_metadata ){
	meta_base = params.input.metadata.get("base", "")
	if(meta_base != "" && !meta_base.endsWith("/") ){
		meta_base = "${meta_base}/"
	}

	Channel.from(params.input.metadata.files.collect{k,v ->
	 	                  [k, file("${meta_base}${v}") ] })
										     .into{METADATA_PRE; METADATA_PRE_VERBOSE }

  BED_FOR_RIBO
                   .join(METADATA_PRE, remainder: true)
                   .into{METADATA_RIBO; METADATA_VERBOSE}
}
else {
  METADATA_PRE = Channel.from([null])
  METADATA_PRE_VERBOSE = Channel.from([null])

  BED_FOR_RIBO
  .map{ sample, bed -> [sample, bed, null] }
  .into{ METADATA_RIBO; METADATA_VERBOSE }
}


if (params.input.get("root_meta", false)){
  ROOT_META = Channel.from([file(params.input.root_meta)])
}
else {
  ROOT_META = Channel.from([null])
}

// METADATA CHANNELS
///////////////////////////////////////////////////////////////////////////////

// FOR DEBUGGING
// QUICK WAY TO CANCEL RIBO CREATION
//METADATA_RIBO = Channel.empty()

////////////////////////////////////////////////////////////////////////////////
/* CREATE RIBO FILES */

do_create_ribo = params.get("do_ribo_creation", true)

if(do_create_ribo){

Channel.from( file(params.input.reference.transcript_lengths) ).
into{T_LENGTHS_FOR_RIBO; T_LENGTHS_FOR_IND_RIBO}

Channel.from( file(params.input.reference.regions) ).
into{ANNOTATION_FOR_RIBO; ANNOTATION_FOR_IND_RIBO}

if(params.ribo.coverage){
	coverage_argument = ""
}
else{
  coverage_argument = "--nocoverage"
}




process create_ribo{

	publishDir get_publishdir("ribo") + "/experiments", mode:'copy'
  storeDir get_storedir("ribo") + "/experiments"

	input:
	set val(sample), file(bed_file), file(meta_file) from METADATA_RIBO
	file(transcript_length_file) from T_LENGTHS_FOR_RIBO.first()
	file(annotation_file) from ANNOTATION_FOR_RIBO.first()
  file(root_meta_file) from ROOT_META.first()

	output:
	set val(sample), file("${sample}.ribo") into RIBO_MAIN

  script:
  if (meta_file != null){
    sample_meta_argument = "--expmeta ${meta_file}"
  }
  else {
    sample_meta_argument = ""
  }

  if(root_meta_file == null){
    root_meta_argument = ""
  }
  else{
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


RIBO_MAIN.into{RIBO_FOR_RNASEQ; RIBO_AFTER_CREATION}

} // end if(do_create_ribo)

// CREATE RIBO FILES
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
/* Post Genome */

do_post_genome = params.input.reference.get("post_genome", false)

if(do_align_genome && do_post_genome ){

  POST_GENOME_INDEX = Channel.from([[
                  params.input.reference.post_genome
                  .split('/')[-1]
                  .replaceAll('\\*$', "")
                  .replaceAll('\\.$', ""),
               file(params.input.reference.post_genome),
              ]])



process post_genome_alignment{

	storeDir get_storedir("post_genome_alignment") + "/" + params.output.individual_lane_directory

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



POST_GENOME_ALIGNMENT_ALIGNED.into{ POST_GENOME_ALIGNMENT_ALIGNED_FASTQ_READ_LENGTH;
                                    POST_GENOME_ALIGNMENT_ALIGNED_MERGE;
                                    POST_GENOME_ALIGNMENT_ALIGNED_FASTQ_FASTQC }

POST_GENOME_ALIGNMENT_UNALIGNED.into{ POST_GENOME_ALIGNMENT_UNALIGNED_FASTQ_READ_LENGTH;
                                      POST_GENOME_ALIGNMENT_UNALIGNED_MERGE;
                                      POST_GENOME_ALIGNMENT_UNALIGNED_FASTQ_FASTQC}

// POST_GENOME ALIGNMENT
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* MERGE POST_GENOME ALIGNMENT */
POST_GENOME_ALIGNMENT_LOG.into{ POST_GENOME_ALIGNMENT_LOG_MERGE; POST_GENOME_ALIGNMENT_LOG_TABLE }


POST_GENOME_ALIGNMENT_BAM.map{sample, index, bam -> [sample, bam]}.groupTuple()
    .set{ POST_GENOME_ALIGNMENT_GROUPED_BAM }

POST_GENOME_ALIGNMENT_ALIGNED_MERGE.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
    .set{ POST_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ }

POST_GENOME_ALIGNMENT_UNALIGNED_MERGE.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
    .set{ POST_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ }
POST_GENOME_ALIGNMENT_LOG_MERGE.map{sample, index, log -> [sample, log]}.groupTuple()
    .set{ POST_GENOME_ALIGNMENT_GROUPED_LOG }


POST_GENOME_ALIGNMENT_GROUPED_BAM.join( POST_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ )
                            .join(POST_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ)
                            .join(POST_GENOME_ALIGNMENT_GROUPED_LOG)
                            .set{ POST_GENOME_ALIGNMENT_GROUPED_JOINT }


process merge_post_genome_alignment{

	storeDir get_storedir("post_genome_alignment") + "/" + params.output.merged_lane_directory

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
   .map{ sample, index, stats_file -> stats_file }
   .toSortedList().set{POST_GENOME_ALIGNMENT_CSV_INDIVIDUAL_LIST}

POST_GENOME_ALIGNMENT_MERGED_CSV
   .map{ sample, stats_file -> stats_file }
   .toSortedList().set{POST_GENOME_ALIGNMENT_CSV_MERGED_LIST}

process combine_individual_postgenome_stats{
  storeDir get_storedir("post_genome_alignment") + "/logs"

  input:
  file(stats_input_files) from POST_GENOME_ALIGNMENT_CSV_INDIVIDUAL_LIST
  file(stats_input_files_merged) from POST_GENOME_ALIGNMENT_CSV_MERGED_LIST

  output:
  file("postgenome_individual_stats.csv") \
        into POST_GENOME_ALIGNMENT_CSV_INDIVIDUAL_COMBINED
  file("postgenome_merged_stats.csv") \
        into POST_GENOME_ALIGNMENT_CSV_MERGED_COMBINED

  """
  rfc merge overall-stats -o postgenome_individual_stats.csv ${stats_input_files} ; \
  rfc merge overall-stats -o postgenome_merged_stats.csv ${stats_input_files_merged}
  """

}

// MERGE POST GENOME ALIGNMENT
////////////////////////////////////////////////////////////////////////////////

} // if( params.input.reference.get("post_genome", false) )

// Post Genome
////////////////////////////////////////////////////////////////////////////////


// TODO: Rethink genome stats approach - temporarily disabled
/*
if(do_align_genome){

  // Append genome stats to transcriptome stats
  process append_genome_stats{
    storeDir get_storedir("stats")

    executor 'local'

    input:
    file(genome_alignment_individual) from GENOME_ALIGNMENT_CSV_INDIVIDUAL_COMBINED
    file(genome_alignment_merged)     from GENOME_ALIGNMENT_CSV_MERGED_COMBINED
    file(individual_alignment_stats)  from COMBINED_INDIVIDUAL_ALIGNMENT_STATS
    file(merged_alignment_stats)      from COMBINED_MERGED_ALIGNMENT_STATS

    output:
    file("individual_stats_with_genome.csv") \
        into COMBINED_INDIVIDUAL_ALIGNMENT_STATS_WITH_GENOME
    file("merged_alignment_stats_with_genome.csv") \
        into COMBINED_MERGED_ALIGNMENT_STATS_WITH_GENOME

    """
    rfc merge concat-csv -o individual_stats_with_genome.csv  \
         ${individual_alignment_stats} ${genome_alignment_individual} && \
    rfc merge concat-csv -o merged_alignment_stats_with_genome.csv \
         ${merged_alignment_stats} ${genome_alignment_merged}
    """
  }

  // Copy comprehensive merged stats to log directory
  process copy_comprehensive_stats_to_log{
    storeDir get_storedir("log") + "/" + params.output.merged_lane_directory

    executor 'local'

    input:
    file(comprehensive_merged_stats) from COMBINED_MERGED_ALIGNMENT_STATS_WITH_GENOME

    output:
    file("*.merged.alignment_stats.csv") into COMPREHENSIVE_MERGED_STATS_IN_LOG

    """
    # Extract sample name from the comprehensive stats file and copy with proper naming
    sample_name=\$(head -1 ${comprehensive_merged_stats} | cut -d',' -f2)
    cp ${comprehensive_merged_stats} \${sample_name}.merged.alignment_stats.csv
    """
  }

  COMBINED_INDIVIDUAL_ALIGNMENT_STATS_WITH_GENOME.set{FINAL_INDIVIDUAL_STATS}
  COMBINED_MERGED_ALIGNMENT_STATS_WITH_GENOME.set{FINAL_MERGED_STATS}
*/

if(do_align_genome && do_align_transcriptome){
  // When both are enabled, use transcriptome stats as final (legacy behavior)
  COMBINED_INDIVIDUAL_ALIGNMENT_STATS.set{FINAL_INDIVIDUAL_STATS}
  COMBINED_MERGED_ALIGNMENT_STATS.set{FINAL_MERGED_STATS}
}
else if(do_align_genome){
  // When only genome is enabled, use genome stats as final
  COMBINED_INDIVIDUAL_GENOME_ALIGNMENT_STATS.set{FINAL_INDIVIDUAL_STATS}
  COMBINED_MERGED_GENOME_ALIGNMENT_STATS.set{FINAL_MERGED_STATS}
}
else if(do_align_transcriptome){
  // When only transcriptome is enabled, use transcriptome stats as final
  COMBINED_INDIVIDUAL_ALIGNMENT_STATS.set{FINAL_INDIVIDUAL_STATS}
  COMBINED_MERGED_ALIGNMENT_STATS.set{FINAL_MERGED_STATS}
}
else{
  // Neither enabled - empty channels
  FINAL_INDIVIDUAL_STATS = Channel.empty()
  FINAL_MERGED_STATS = Channel.empty()
}

// Final stats channels are set above in the genome/transcriptome conditional blocks

// Use final stats directly for publishing
FINAL_INDIVIDUAL_STATS.set{FINAL_INDIVIDUAL_STATS_FOR_PUBLISH}
FINAL_MERGED_STATS.set{FINAL_MERGED_STATS_FOR_PUBLISH}

// SUMMARY STATS - Combine transcriptome and genome when both are available
///////////////////////////////////////////////////////////////////////////////

// Summary stats processes temporarily removed until genome CSV generation is reimplemented
// if(do_align_transcriptome && do_align_genome){
//   process create_summary_individual_stats{
//     storeDir get_storedir("stats") + "/summary/individual"
//
//     input:
//     file(transcriptome_stats) from COMBINED_INDIVIDUAL_ALIGNMENT_STATS
//     file(genome_stats) from GENOME_ALIGNMENT_CSV_INDIVIDUAL_COMBINED
//
//     output:
//     file("summary_individual_stats.csv") \
//           into SUMMARY_INDIVIDUAL_STATS_COMBINED
//
//     """
//     rfc merge concat-csv -o summary_individual_stats.csv \
//          ${transcriptome_stats} ${genome_stats}
//     """
//
//   }
//
//   process create_summary_merged_stats{
//     storeDir get_storedir("stats") + "/summary/merged"
//
//     input:
//     file(transcriptome_stats) from COMBINED_MERGED_ALIGNMENT_STATS
//     file(genome_stats) from GENOME_ALIGNMENT_CSV_MERGED_COMBINED
//
//     output:
//     file("summary_merged_stats.csv") \
//           into SUMMARY_MERGED_STATS_COMBINED
//
//     """
//     rfc merge concat-csv -o summary_merged_stats.csv \
//          ${transcriptome_stats} ${genome_stats}
//     """
//
//   }
// }

process publish_stats{

  publishDir get_publishdir("stats"), mode: "copy"

  executor 'local'

  input:
  file(individual_stats) from FINAL_INDIVIDUAL_STATS_FOR_PUBLISH
  file(merged_stats)     from FINAL_MERGED_STATS_FOR_PUBLISH

  output:
  file("individual_stats.csv") into INDIVIDUAL_STATS_PUBLISHED
  file("stats.csv")            into MERGED_STATS_PUBLISHED

  """
  cp ${individual_stats} individual_stats.csv && \
  cp ${merged_stats} stats.csv
  """
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////                       /* RNA-Seq */                            /////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////// General Function Definitions ////////////////////////////////////////////

String get_rnaseq_storedir(output_type){
    new File( params.output.intermediates.base + "/rnaseq",
              params.output.intermediates.get(output_type, output_type) )
							.getCanonicalPath()
}

String get_rnaseq_publishdir(output_type){
    new File( params.output.output.base + "/rnaseq",
              params.output.output.get(output_type, output_type) )
							.getCanonicalPath()
}

////// General Function Definitions ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



// Both the boolean flag 'do_rnaseq'
// AND actual rnaseq node must be set to perform
// rnaseq data processing steps.
do_rnaseq = params.get("do_rnaseq", false) && \
            params.get("rnaseq", false)

// This outer if clause contains the rest of the RNASEQ
if (do_rnaseq){

rnaseq_fastq_base =  params.rnaseq.get("fastq_base", "")
if(! rnaseq_fastq_base.endsWith("/") && rnaseq_fastq_base != "") {
	rnaseq_fastq_base = "${rnaseq_fastq_base}/"
}

// Group input files into a list of tuples where each item is
// [ sample, fileindex, path_to_fastq_file]

Channel.from(params.rnaseq.fastq.collect{k,v ->
	              v.collect{ z -> [k, v.indexOf(z) + 1,
									               file("${rnaseq_fastq_base}${z}")] }  })
	.flatten().collate(3).into{  RNASEQ_FASTQ;
                               RNASEQ_FASTQ_VERBOSE;
                               RNASEQ_FASTQ_FASTQC;
                               RNASEQ_FASTQ_CLIP;
                               RNASEQ_FASTQ_EXISTENCE}

if(params.do_check_file_existence){
  // Make Sure Fastq Files Exist
  RNASEQ_FASTQ_EXISTENCE
  .map{ sample, index, this_file -> file_exists(this_file) }
}

process rnaseq_raw_fastqc{

  	publishDir get_rnaseq_publishdir("fastqc"), mode: 'copy'

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

process rnaseq_clip{
  storeDir get_rnaseq_storedir("clip")

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

RNASEQ_FILTER_INDEX = Channel.from([[
             params.input.reference.filter
                .split('/')[-1]
                .replaceAll('\\*$', "")
                .replaceAll('\\.$', ""),
             file(params.input.reference.filter),
            ]])

process rnaseq_filter{

	storeDir get_rnaseq_storedir("filter")

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

RNASEQ_FILTER_UNALIGNED.into{RNASEQ_FILTER_UNALIGNED_FASTQ_READ_LENGTH;
                             RNASEQ_FILTER_UNALIGNED_FASTQ_FASTQC;
                             RNASEQ_FILTER_UNALIGNED_TRANSCRIPTOME;
                             RNASEQ_FILTER_UNALIGNED_GENOME}


rnaseq_bt2_arguments = params.rnaseq.get("bt2_argumments", "")

RNASEQ_TRANSCRIPTOME_INDEX = Channel.from([[
            params.input.reference.transcriptome
               .split('/')[-1]
               .replaceAll('\\*$', "")
               .replaceAll('\\.$', ""),
            file(params.input.reference.transcriptome),
           ]])


process rnaseq_transcriptome_alignment{

   storeDir get_rnaseq_storedir("transcriptome_alignment") + "/" +\
                                 params.output.individual_lane_directory

   input:
   set val(sample), val(index), file(fastq) \
                from RNASEQ_FILTER_UNALIGNED_TRANSCRIPTOME
	 set val(transcriptome_reference), file(transcriptome_Reference_files) \
		            from RNASEQ_TRANSCRIPTOME_INDEX.first()

   output:
   set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.bam") \
       into RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAM_PRE
   set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.bam.bai") \
       into RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAI
   set val(sample), val(index), file("${sample}.${index}.aligned.transcriptome_alignment.fastq.gz") \
       into RNASEQ_TRANSCRIPTOME_ALIGNMENT_ALIGNED
   set val(sample), val(index), file("${sample}.${index}.unaligned.transcriptome_alignment.fastq.gz") \
       into RNASEQ_TRANSCRIPTOME_ALIGNMENT_UNALIGNED
   set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.log") \
       into RNASEQ_TRANSCRIPTOME_ALIGNMENT_LOG
   set val(sample), val(index), file("${sample}.${index}.transcriptome_alignment.stats") \
       into RNASEQ_TRANSCRIPTOME_ALIGNMENT_STATS

   """
   bowtie2 ${rnaseq_bt2_arguments} \
           -x ${transcriptome_reference} -q ${fastq} \
           --threads ${task.cpus} \
           --al-gz ${sample}.${index}.aligned.transcriptome_alignment.fastq.gz \
           --un-gz ${sample}.${index}.unaligned.transcriptome_alignment.fastq.gz \
						           2> ${sample}.${index}.transcriptome_alignment.log \
           | samtools view -bS - \
           | samtools sort -@ ${task.cpus} -o ${sample}.${index}.transcriptome_alignment.bam \
           && samtools index -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.bam \
           && samtools idxstats -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.bam  > \
              ${sample}.${index}.transcriptome_alignment.stats
   """
}



RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAM_PRE
.into{ RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAM;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAM_MERGE;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAM_FOR_QUALITY}

process rnaseq_quality_filter{

	storeDir get_rnaseq_storedir("quality_filter")

	input:
	set val(sample), val(index), file(bam) \
        from RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAM_FOR_QUALITY

	output:
	set val(sample), val(index),
        file("${sample}.${index}.transcriptome_alignment.qpass.bam") \
        into RNASEQ_TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_PRE
    set val(sample), val(index),
        file("${sample}.${index}.transcriptome_alignment.qpass.bam.bai") \
        into RNASEQ_TRANSCRIPTOME_ALIGNMENT_QPASS_BAI
    set val(sample), val(index),
        file("${sample}.${index}.qpass.count") \
        into RNASEQ_TRANSCRIPTOME_QPASS_COUNTS
    set val(sample), val(index),
        file("${sample}.${index}.transcriptome_alignment.qpass.stats") \
        into RNASEQ_TRANSCRIPTOME_ALIGNMENT_QPASS_STATS

	"""
	samtools view -b -q ${params.mapping_quality_cutoff} ${bam}\
	| samtools sort -@ ${task.cpus} -o ${sample}.${index}.transcriptome_alignment.qpass.bam \
	&& samtools view -b -c ${sample}.${index}.transcriptome_alignment.qpass.bam > ${sample}.${index}.qpass.count \
	&& samtools index -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.qpass.bam \
	&& samtools idxstats -@ ${task.cpus} ${sample}.${index}.transcriptome_alignment.qpass.bam  > \
               ${sample}.${index}.transcriptome_alignment.qpass.stats
	"""
}

RNASEQ_TRANSCRIPTOME_ALIGNMENT_QPASS_BAM_PRE
.into{ RNASEQ_QPASS_BAM_READ_LENGTH;
	     RNASEQ_TRANSCRIPTOME_ALIGNMENT_QPASS_BAM}

// QUALITY FILTER
///////////////////////////////////////////////////////////////////////////////////////

RNASEQ_TRANSCRIPTOME_QPASS_COUNTS
.into{RNASEQ_TRANSCRIPTOME_QPASS_COUNTS_FOR_INDEX;
	    RNASEQ_TRANSCRIPTOME_QPASS_COUNTS_FOR_TABLE}

// We need to copy output channels of transcriptome alignment
// for merging and variaous steps of downstream processing

RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAI
.into{ RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAI_MERGE ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_BAI_REGION_COUNT}

RNASEQ_TRANSCRIPTOME_ALIGNMENT_ALIGNED
.into{ RNASEQ_TRANSCRIPTOME_ALIGNMENT_ALIGNED_MERGE ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_ALIGNED_LENGTH ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_ALIGNED_FASTQC }

RNASEQ_TRANSCRIPTOME_ALIGNMENT_UNALIGNED
.into{ RNASEQ_TRANSCRIPTOME_ALIGNMENT_UNALIGNED_MERGE ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_UNALIGNED_GENOME ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_UNALIGNED_LENGTH ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_UNALIGNED_FASTQC }

RNASEQ_TRANSCRIPTOME_ALIGNMENT_LOG
.into{ RNASEQ_TRANSCRIPTOME_ALIGNMENT_LOG_MERGE ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_LOG_TABLE  }

RNASEQ_TRANSCRIPTOME_ALIGNMENT_STATS
.into{ RNASEQ_TRANSCRIPTOME_ALIGNMENT_STATS_MERGE ;
       RNASEQ_TRANSCRIPTOME_ALIGNMENT_STATS_TABLE  }


process rnaseq_bam_to_bed{

	storeDir get_rnaseq_storedir("bam_to_bed") + "/" + params.output.individual_lane_directory

	input:
	set val(sample), val(index), file(bam) from RNASEQ_TRANSCRIPTOME_ALIGNMENT_QPASS_BAM

	output:
	set val(sample), val(index), file("${sample}.${index}.bed") into RNASEQ_BAM_TO_BED
	set val(sample), val(index), file("${sample}.${index}_nodedup_count.txt") \
	   into RNASEQ_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP

   """
   if [ `samtools view -c ${bam}` -eq 0 ];
   then
      touch ${sample}.${index}.bed
   else
       bamToBed -i ${bam} > ${sample}.${index}.bed
   fi

   wc -l ${sample}.${index}.bed > ${sample}.${index}_nodedup_count.txt
   """
}

 RNASEQ_BAM_TO_BED.into{  RNASEQ_BED_NODEDUP;
                          RNASEQ_BED_FOR_DEDUP;
                          RNASEQ_BED_FOR_INDEX_SEP_PRE }

// RNA-seq deduplication now controlled by universal dedup_method
// Remove separate RNA-seq deduplicate flag - use dedup_method != "none"

process rnaseq_add_sample_index_col_to_bed{

	storeDir get_rnaseq_storedir("bam_to_bed") + "/" + params.output.individual_lane_directory

	input:
   set val(sample), val(index), file(bed) from  RNASEQ_BED_FOR_DEDUP

	output:
	set val(sample), file("${sample}.${index}.with_sample_index.bed")\
	     into  RNASEQ_BED_FOR_DEDUP_INDEX_COL_ADDED

	"""
	awk -v newcol=${sample}.${index} '{print(\$0"\\t"newcol)}' ${bed}\
	   > ${sample}.${index}.with_sample_index.bed
	"""
}

RNASEQ_BED_FOR_DEDUP_INDEX_COL_ADDED.groupTuple()
   .set{  RNASEQ_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED }

process rnaseq_merge_bed{

	storeDir get_rnaseq_storedir("bam_to_bed") + "/" + params.output.merged_lane_directory

	input:
	set val(sample), file(bed_files) from  RNASEQ_BED_FOR_DEDUP_INDEX_COL_ADDED_GROUPED

	output:
	set val(sample), file("${sample}.merged.pre_dedup.bed") \
	    into  RNASEQ_BED_MERGED_PRE_DEDUP


	"""
	cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.merged.pre_dedup.bed
	"""
}

RNASEQ_BED_MERGED_PRE_DEDUP
.into{RNASEQ_BED_FOR_DEDUP_MERGED_PRE_DEDUP;
      RNASEQ_BED_NODEDUP_FOR_RIBO}


process rnaseq_deduplicate{

	storeDir  get_rnaseq_storedir("alignment_ribo") + "/" + params.output.merged_lane_directory

	input:
	set val(sample), file(bed) from RNASEQ_BED_FOR_DEDUP_MERGED_PRE_DEDUP

	output:
	set val(sample), file("${sample}.merged.post_dedup.bed") \
	     into RNASEQ_BED_FOR_DEDUP_MERGED_POST_DEDUP

	when:
	dedup_method != "none"

	"""
	rfc dedup -i ${bed} -o ${sample}.merged.post_dedup.bed
	"""
}

RNASEQ_BED_FOR_DEDUP_MERGED_POST_DEDUP
.into{RNASEQ_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_SEP;
      RNASEQ_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_RIBO}

RNASEQ_BED_FOR_INDEX_SEP_PRE
.map{ sample,index,file -> [sample, index] }
.combine(RNASEQ_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_SEP, by:0)
.set{ RNASEQ_BED_FOR_INDEX_SEP_POST_DEDUP }

process rnaseq_separate_bed_post_dedup{

	storeDir  get_rnaseq_storedir("alignment_ribo") + "/" + params.output.individual_lane_directory

	input:
	set val(sample), val(index), file(bed) from  RNASEQ_BED_FOR_INDEX_SEP_POST_DEDUP

	output:
	set val(sample), val(index), file("${sample}.${index}.post_dedup.bed") \
	   into  RNASEQ_BED_DEDUPLICATED
	set val(sample), val(index), file("${sample}.${index}.count_after_dedup.txt")\
	   into  RNASEQ_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP

	"""
	awk -v this_sample=${sample}.${index} \
	 '{ if(\$7 == this_sample ){print(\$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6)} }' ${bed} > ${sample}.${index}.post_dedup.bed \
	  && wc -l ${sample}.${index}.post_dedup.bed > ${sample}.${index}.count_after_dedup.txt
	"""
}

if(dedup_method != "none"){
  RNASEQ_BED_FOR_DEDUP_MERGED_POST_DEDUP_FOR_RIBO
     .into{RNASEQ_BED_FOR_SEPARATION; RNASEQ_BED_FOR_RIBO_FINAL}
  RNASEQ_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP.set{RNASEQ_INDIVIDUAL_DEDUP_COUNT}
}
else{
  RNASEQ_BED_NODEDUP_FOR_RIBO
      .into{RNASEQ_BED_FOR_SEPARATION; RNASEQ_BED_FOR_RIBO_FINAL}
  RNASEQ_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP.set{RNASEQ_INDIVIDUAL_DEDUP_COUNT}
}


////////////////////////////////////////////////////////////////////////////////


// We need to group the log files by sample name and index
// than flatten that list and group again so that each
// entry can be emmited in groups of 6 for each task


RNASEQ_CLIP_LOG.map{ sample, index, clip_log -> [ [sample, index], clip_log ] }
        .set{RNASEQ_CLIP_LOG_INDEXED}
RNASEQ_FILTER_LOG.map{ sample, index, filter_log -> [ [sample, index], filter_log ] }
          .set{RNASEQ_FILTER_LOG_INDEXED}
RNASEQ_TRANSCRIPTOME_ALIGNMENT_LOG_TABLE
    .map{ sample, index, transcriptome_log -> [ [sample, index], transcriptome_log ] }
    .set{RNASEQ_TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED}
RNASEQ_TRANSCRIPTOME_QPASS_COUNTS_FOR_INDEX
    .map{ sample, index, qpass_count -> [ [sample, index], qpass_count ] }
    .set{RNASEQ_TRANSCRIPTOME_QPASS_COUNTS_INDEXED}
RNASEQ_INDIVIDUAL_DEDUP_COUNT
     .map{ sample, index, dedup_count -> [ [sample, index], dedup_count ] }
     .set{ RNASEQ_INDIVIDUAL_DEDUP_COUNT_INDEXED }

RNASEQ_CLIP_LOG_INDEXED.join(RNASEQ_FILTER_LOG_INDEXED)
                .join(RNASEQ_TRANSCRIPTOME_ALIGNMENT_LOG_TABLE_INDEXED)
                .join(RNASEQ_TRANSCRIPTOME_QPASS_COUNTS_INDEXED)
                .join(RNASEQ_INDIVIDUAL_DEDUP_COUNT_INDEXED)
                .flatten()
                .collate(7)
                .set{ RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_INPUT }


process rnaseq_individual_alignment_stats{

	//Compiles statistics coming from the individual steps:
	//cutadapt, filter, transcriptome and genome alignment,
	//quality filtering and deduplication


	 executor 'local'

   storeDir get_rnaseq_storedir("stats")

   input:
   set val(sample), val(index), file(clip_log), file(filter_log),\
       file(transcriptome_log), file(qpass_count),\
       file(dedup_count)\
       from RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_INPUT

   output:
   set val(sample), val(index), file("${sample}.${index}.rnaseq_overall_alignment.csv") \
      into RNASEQ_INDIVIDUAL_ALIGNMENT_STATS

   """
   rfc compile-step-stats \
	   -n ${sample}.${index} \
     -c ${clip_log} \
     -f ${filter_log} \
     -t ${transcriptome_log} \
     -q ${qpass_count} \
     -d ${dedup_count} \
     -o ${sample}.${index}.rnaseq_overall_alignment.csv
   """

}

// INDIVIDUAL ALIGNMENT STATS
///////////////////////////////////////////////////////////////////////////////////////

RNASEQ_INDIVIDUAL_ALIGNMENT_STATS
    .into{ RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION;
           RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING}

///////////////////////////////////////////////////////////////////////////////////////
/* COMBINE INDIVIDUAL ALIGNMENT STATS */

RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_FOR_COLLECTION
    .map{ sample, index, stats_file -> stats_file }
    .toSortedList().set{RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED}

process rnaseq_combine_individual_alignment_stats{

  executor 'local'

	publishDir get_rnaseq_publishdir("stats"), mode: 'copy'

	input:
	file(stat_table) from RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_COLLECTED

	output:
	file("rnaseq_individual_stats.csv") \
	      into RNASEQ_COMBINED_INDIVIDUAL_ALIGNMENT_STATS

	"""
	  rfc merge overall-stats \
	   -o raw_combined_individual_aln_stats.csv \
	      ${stat_table} && \
    rfc stats-percentage \
	  -i raw_combined_individual_aln_stats.csv \
	  -o rnaseq_individual_stats.csv
	"""
}


// COMBINE INDIVIDUAL ALIGNMENT STATS
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
/* SUM INDIVIDUAL ALIGNMENT STATS */

/*
For each sample, sums up the stats coming from individual lanes
*/

RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_FOR_GOUPING
    .map{ sample, index, file -> [ sample, file ] }
    .groupTuple()
    .into{ RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_GROUPED ;
           RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_GROUPED_VERBOSE }

process rnaseq_sum_individual_alignment_stats{

  executor 'local'

	storeDir get_rnaseq_storedir( "log/" + params.output.merged_lane_directory )

	input:
	set val(sample), file(stat_files) from RNASEQ_INDIVIDUAL_ALIGNMENT_STATS_GROUPED

	output:
	set val(sample), file("${sample}.rnaseq.merged.alignment_stats.csv")\
	   into RNASEQ_MERGED_ALIGNMENT_STATS

	"""
	rfc sum-stats -n ${sample}\
	  -o ${sample}.rnaseq.merged.alignment_stats.csv ${stat_files}
	"""
}

// SUM INDIVIDUAL ALIGNMENT STATS
////////////////////////////////////////////////////////////////////////////////

RNASEQ_MERGED_ALIGNMENT_STATS
.map{ sample, stats_file -> stats_file }
.toSortedList()
.set{ RNASEQ_MERGED_ALIGNMENT_STATS_COLLECTED }

////////////////////////////////////////////////////////////////////////////////
/* COMBINE MERGED ALIGNMENT STATS */

process rnaseq_combine_merged_alignment_stats{

	publishDir get_rnaseq_publishdir("stats"), mode: 'copy'

	executor 'local'

	input:
	file(stat_files) from RNASEQ_MERGED_ALIGNMENT_STATS_COLLECTED

	output:
	file("rnaseq_stats.csv") into RNASEQ_COMBINED_MERGED_ALIGNMENT_STATS

	"""
	rfc merge overall-stats \
	    -o raw_combined_merged_aln_stats.csv \
	    ${stat_files} && \
	rfc stats-percentage \
	  -i raw_combined_merged_aln_stats.csv \
	  -o rnaseq_stats.csv
	"""
}

// COMBINE MERGED ALIGNMENT STATS
////////////////////////////////////////////////////////////////////////////////

RNASEQ_FOR_RIBOPY   = Channel.create()
RIBO_FOR_RNASEQ_EXCLUDED = Channel.create()

/*
Separate the ribo files which have rnaseq data and which dont
*/
RIBO_FOR_RNASEQ
.join(RNASEQ_BED_FOR_RIBO_FINAL, remainder: true)
.choice( RNASEQ_FOR_RIBOPY, RIBO_FOR_RNASEQ_EXCLUDED )
{it[2] != null ? 0 : 1}

RIBO_FOR_RNASEQ_EXCLUDED.map{ sample, ribo, bed_null -> [sample, ribo]}
.set{ RIBO_FOR_RNASEQ_EXCLUDED_FOR_MERGE }

process put_rnaseq_into_ribo{
  publishDir get_publishdir("ribo") + "/experiments", mode: 'copy'

  input:
  set val(sample), file(ribo), file(rnaseq) from RNASEQ_FOR_RIBOPY

  output:
  set val(sample), file(ribo) into RIBO_WITH_RNASEQ_PRE

  """
  ribopy rnaseq set -n ${sample} -a ${rnaseq} -f bed --force ${ribo}
  """

}

// For the downsrtream "merge_ribo" process,
// we need to combine the ribos with and without rnaseq data.
RIBO_WITH_RNASEQ_PRE.concat( RIBO_FOR_RNASEQ_EXCLUDED_FOR_MERGE )
.into{RIBO_WITH_RNASEQ; RIBO_WITH_RNASEQ_VERBOSE}

///////////////////////////////////////////////////////////////////////////////
//          R N A - S E Q   G E N O M E   A L I G N M E N T
///////////////////////////////////////////////////////////////////////////////

// RNA-seq genome alignment (conditional on genome alignment being enabled)
if(do_align_genome) {

    // Set up RNA-seq genome input channel (similar to ribo-seq)
    if(do_align_transcriptome) {
        RNASEQ_GENOME_INPUT_CHANNEL = RNASEQ_TRANSCRIPTOME_ALIGNMENT_UNALIGNED_GENOME
    } else {
        RNASEQ_GENOME_INPUT_CHANNEL = RNASEQ_FILTER_UNALIGNED_GENOME
    }

    process rnaseq_genome_alignment{

        storeDir get_rnaseq_storedir("genome_alignment") + "/" + params.output.individual_lane_directory

        input:
        set val(sample), val(index), file(fastq) from RNASEQ_GENOME_INPUT_CHANNEL
        set val(genome_base), file(genome_files) from GENOME_INDEX.first()

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
        hisat2 ${params.rnaseq.get("hisat2_arguments", params.alignment_arguments.genome)} \\
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
               && python3 ${workflow.projectDir}/hisat2-log-to-csv.py \\
                      -l ${sample}.${index}.rnaseq_genome_alignment.log \\
                      -n ${sample} -p rnaseq_genome \\
                      -o ${sample}.${index}.rnaseq_genome_alignment.csv
        """
    }

    // Group RNA-seq genome alignment outputs for merging
    RNASEQ_GENOME_ALIGNMENT_BAM.map{sample, index, bam -> [sample, bam]}.groupTuple()
        .set{ RNASEQ_GENOME_ALIGNMENT_GROUPED_BAM }

    RNASEQ_GENOME_ALIGNMENT_ALIGNED.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
        .set{ RNASEQ_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ }

    RNASEQ_GENOME_ALIGNMENT_UNALIGNED.map{sample, index, fastq -> [sample, fastq]}.groupTuple()
        .set{ RNASEQ_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ }

    RNASEQ_GENOME_ALIGNMENT_LOG.map{sample, index, log -> [sample, log]}.groupTuple()
        .set{ RNASEQ_GENOME_ALIGNMENT_GROUPED_LOG }

    // Join all grouped channels for merging
    RNASEQ_GENOME_ALIGNMENT_GROUPED_BAM.join( RNASEQ_GENOME_ALIGNMENT_GROUPED_ALIGNED_FASTQ )
                                .join(RNASEQ_GENOME_ALIGNMENT_GROUPED_UNALIGNED_FASTQ)
                                .join(RNASEQ_GENOME_ALIGNMENT_GROUPED_LOG)
                                .set{ RNASEQ_GENOME_ALIGNMENT_GROUPED_JOINT }

    process rnaseq_merge_genome_alignment{

        storeDir get_rnaseq_storedir("genome_alignment") + "/" + params.output.merged_lane_directory

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
        cp ${alignment_log} ${sample}.rnaseq_genome.log && \\
        python3 ${workflow.projectDir}/hisat2-log-to-csv.py \\
                           -l ${sample}.rnaseq_genome.log -n ${sample} -p rnaseq_genome \\
                           -o ${sample}.rnaseq_genome.csv
        """
    }

    ///////////////////////////////////////////////////////////////////////////////
    //          R N A - S E Q   G E N O M E   D E D U P L I C A T I O N
    ///////////////////////////////////////////////////////////////////////////////

    // Split RNA-seq genome BAM channel for deduplication
    RNASEQ_GENOME_ALIGNMENT_MERGED_BAM.into{
        RNASEQ_GENOME_BAM_FOR_DEDUP;
        RNASEQ_GENOME_BAM_NODEDUP
    }

    process rnaseq_genome_bam_to_bed{

        storeDir get_rnaseq_storedir("bam_to_bed") + "/" + params.output.merged_lane_directory

        input:
        set val(sample), file(bam) from RNASEQ_GENOME_BAM_FOR_DEDUP

        output:
        set val(sample), file("${sample}.rnaseq_genome.bed") into RNASEQ_GENOME_BED_FOR_DEDUP_PRE
        set val(sample), file("${sample}.rnaseq_genome_nodedup_count.txt") \
           into RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP

        when:
        dedup_method != "none"

        """
        if [ `samtools view -c ${bam}` -eq 0 ];
        then
           touch ${sample}.rnaseq_genome.bed
        else
            bamToBed -i ${bam} > ${sample}.rnaseq_genome.bed
        fi

        wc -l ${sample}.rnaseq_genome.bed > ${sample}.rnaseq_genome_nodedup_count.txt
        """
    }

    // RNA-seq genome deduplication (position-based only)
    process rnaseq_genome_deduplicate{

        storeDir get_rnaseq_storedir("alignment_ribo") + "/" + params.output.merged_lane_directory

        input:
        set val(sample), file(bed) from RNASEQ_GENOME_BED_FOR_DEDUP_PRE

        output:
        set val(sample), file("${sample}.rnaseq_genome.post_dedup.bed") \
             into RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP

        when:
        dedup_method != "none"

        """
        # RNA-seq uses position-based deduplication only (no UMI support)
        rfc dedup -i ${bed} -o ${sample}.rnaseq_genome.post_dedup.bed
        """
    }

    // Create count file for deduplicated RNA-seq genome reads
    process rnaseq_genome_count_deduplicated_reads{
        storeDir get_rnaseq_storedir("alignment_ribo") + "/" + params.output.merged_lane_directory

        input:
        set val(sample), file(bed) from RNASEQ_GENOME_BED_FOR_DEDUP_MERGED_POST_DEDUP

        output:
        set val(sample), file("${sample}.rnaseq_genome_dedup_count.txt") \
            into RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP

        when:
        dedup_method != "none"

        """
        wc -l ${bed} > ${sample}.rnaseq_genome_dedup_count.txt
        """
    }

    // Select the appropriate count channel based on deduplication method
    if(dedup_method != "none"){
      RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITH_DEDUP.set{RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT}
    }
    else{
      RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT_WITHOUT_DEDUP.set{RNASEQ_GENOME_INDIVIDUAL_DEDUP_COUNT}
    }

    ///////////////////////////////////////////////////////////////////////////////
    //          E N D   R N A - S E Q   G E N O M E   D E D U P L I C A T I O N
    ///////////////////////////////////////////////////////////////////////////////

} // end RNA-seq genome alignment conditional

} // if (do_rnaseq)
// RNA-Seq
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Merge Ribos*/

if(do_create_ribo){
  if(do_rnaseq){
    RIBO_WITH_RNASEQ.into{RIBO_FOR_MERGE_PRE; RIBO_FOR_COUNT}
  }
  else{
    RIBO_AFTER_CREATION.into{RIBO_FOR_MERGE_PRE; RIBO_FOR_COUNT}
  }
  RIBO_FOR_MERGE_PRE.map{ sample, ribo -> [ribo]}.flatten().collect()
                    .set{RIBO_FOR_MERGE}

  process merge_ribos{

    publishDir get_publishdir("ribo"), mode:'copy'

    input:
    file(sample_ribo) from RIBO_FOR_MERGE
    val(ribo_count) from RIBO_FOR_COUNT.count()

    output:
    file("all.ribo") into ALL_RIBO

    script:
    if(ribo_count > 1){
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
