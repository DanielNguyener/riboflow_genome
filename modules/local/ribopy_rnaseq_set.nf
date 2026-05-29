// Merge the RNA-seq transcriptome BED into a ribo-seq .ribo file via
// `ribopy rnaseq set`. Operates on a copy so the original .ribo from RIBOPY_CREATE
// is preserved in the work dir; the publishDir saveAs strips `.updated` so the
// published .ribo is overwritten with the RNA-seq-augmented version.
process RIBOPY_RNASEQ_SET {
    tag "${meta.id}"

    input:
    tuple val(meta), path(ribo), path(rnaseq_bed)

    output:
    tuple val(meta), path("${meta.id}.updated.ribo"), emit: ribo

    script:
    """
    cp ${ribo} ${meta.id}.updated.ribo
    cut -f1-6 ${rnaseq_bed} > rnaseq_input.bed
    ribopy rnaseq set -n ${meta.id} -a rnaseq_input.bed -f bed --force ${meta.id}.updated.ribo
    """

    stub:
    """
    touch ${meta.id}.updated.ribo
    """
}
