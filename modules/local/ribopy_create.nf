// Creates a per-sample .ribo file from a merged post-dedup BED.
// BED is stripped to 6 columns before passing to ribopy: the position-dedup
// path appends a 7th sample.lane column for dedup; ribopy only wants cols 1-6.
// All tuning flags come from params.ribo.* (set in test.config / example YAMLs).
process RIBOPY_CREATE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(dedup_bed), path(expmeta_file)
    path regions_bed
    path lengths_tsv

    output:
    tuple val(meta), path("${meta.id}.ribo"), emit: ribo

    script:
    def ref_name     = params.ribo?.ref_name         ?: 'appris_human'
    def radius       = params.ribo?.metagene_radius  ?: 50
    def left_span    = params.ribo?.left_span        ?: 35
    def right_span   = params.ribo?.right_span       ?: 15
    def len_min      = params.ribo?.read_length?.min ?: 15
    def len_max      = params.ribo?.read_length?.max ?: 35
    def nocov_flag   = (params.ribo?.coverage != false) ? '' : '--nocoverage'
    // Per-sample expmeta takes precedence; fall back to global params.ribo.expmeta.
    // Global paths are resolved to absolute so ribopy can find them from the work dir.
    def expmeta_arg  = expmeta_file
        ? "--expmeta ${expmeta_file}"
        : (params.ribo?.expmeta ? "--expmeta ${new File(params.ribo.expmeta.toString()).absolutePath}" : '')
    def ribometa_arg = params.ribo?.ribometa
        ? "--ribometa ${new File(params.ribo.ribometa.toString()).absolutePath}"
        : ''
    """
    cut -f1-6 ${dedup_bed} > ribo_input.bed
    ribopy create \\
        --name ${meta.id} \\
        --alignmentfile ribo_input.bed \\
        --reference ${ref_name} \\
        --lengths ${lengths_tsv} \\
        --annotation ${regions_bed} \\
        --metageneradius ${radius} \\
        -l ${left_span} -r ${right_span} \\
        --lengthmin ${len_min} --lengthmax ${len_max} \\
        ${nocov_flag} \\
        ${expmeta_arg} \\
        ${ribometa_arg} \\
        -n ${task.cpus} \\
        ${meta.id}.ribo
    """

    stub:
    """
    touch ${meta.id}.ribo
    """
}
