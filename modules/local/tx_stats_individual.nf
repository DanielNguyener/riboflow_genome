// Per-lane transcriptome alignment stats row. Parallel to stats_individual.nf
// (genome). Parses the bowtie2 transcriptome log + clip/filter logs + qpass/dedup
// counts into a raw-count CSV. Includes clip/filter rows so transcriptome-only
// runs produce a self-contained stats file.
process TX_STATS_INDIVIDUAL {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta),
          path(tx_log), path(clip_log), path(filter_log),
          path(qpass_total), path(qpass_primary), path(qpass_secondary),
          path('dedup.total.count'), path('dedup.primary.count'), path('dedup.secondary.count')

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_individual.csv"), emit: csv

    script:
    def prefix = "${meta.id}.${meta.lane}"
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

# filter log (bowtie2 rRNA/tRNA)
fl = [l for l in open('${filter_log}') if l.strip() and not l[0].isalpha()]
filtered_out = int(fl[3].split()[0]) + int(fl[4].split()[0])
filter_kept  = int(fl[2].split()[0])

# transcriptome bowtie2 alignment log
tl = [l for l in open('${tx_log}') if l.strip() and not l[0].isalpha()]
tx_unaligned    = int(tl[2].split()[0])
tx_aligned_once = int(tl[3].split()[0])
tx_aligned_many = int(tl[4].split()[0])
tx_primary      = tx_aligned_once + tx_aligned_many

def read_int(p):
    with open(p) as fh:
        return int(fh.read().strip().split()[0])

qpass_total_v = read_int('${qpass_total}')
dedup_total_v = read_int('dedup.total.count')

rows = [
    ('total_reads',                         total_reads),
    ('clipped_reads',                        clipped_reads),
    ('filtered_out',                         filtered_out),
    ('filter_kept',                          filter_kept),
    ('transcriptome_aligned_once',           tx_aligned_once),
    ('transcriptome_aligned_many',           tx_aligned_many),
    ('transcriptome_total_aligned',          tx_primary),
    ('transcriptome_unaligned',              tx_unaligned),
    ('transcriptome_qpass_aligned_reads',    qpass_total_v),
    ('transcriptome_after_dedup',            dedup_total_v),
]
with open('${prefix}.transcriptome_individual.csv', 'w') as fh:
    fh.write(',${prefix}\\n')
    for k, v in rows:
        fh.write(f'{k},{v}\\n')
PYEOF
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    printf ',${prefix}\\ntotal_reads,0\\n' > ${prefix}.transcriptome_individual.csv
    """
}
