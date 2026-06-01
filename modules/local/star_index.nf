// STAR genome index generation (genomeGenerate mode).
// Triggered when the user provides genome_fasta + gtf instead of a pre-built
// index directory. Output directory is piped directly into STAR_ALIGN via
// ch_genome_index = STAR_INDEX.out.index.first() in workflows/riboflow.nf.
process STAR_INDEX {
    tag "star_index"

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    def overhang = task.ext.sjdb_overhang ?: (params.star?.sjdb_overhang ?: 100)
    def extra    = task.ext.args          ?: (params.star?.index_args ?: '')
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${overhang} \\
        ${extra}
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/SA star_index/SAindex star_index/Genome star_index/chrNameLength.txt
    """
}
