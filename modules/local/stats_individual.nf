// Per-lane alignment-stats row. Ports `individual_genome_alignment_stats`
// (RiboFlow.groovy:1528-1617): parses the cutadapt + bowtie2 logs, the STAR log
// via `rfc parse-star-log`, and the qpass/dedup count files into a raw-count CSV.
// dedup_* counts are staged under fixed names so they never collide with the
// qpass_* files when dedup_method=='none' aliases qpass counts into the dedup slot.
process STATS_INDIVIDUAL {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta),
          path(clip_log), path(filter_log), path(genome_log), path(genome_secondary_count),
          path(qpass_total), path(qpass_primary), path(qpass_secondary),
          path('dedup.total.count'), path('dedup.primary.count'), path('dedup.secondary.count')

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.genome_individual.csv"), emit: csv

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
with open('${prefix}.genome_individual.csv', 'w') as fh:
    fh.write(',${prefix}\\n')
    for k, v in rows:
        fh.write(f'{k},{v}\\n')
PYEOF
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    printf ',${prefix}\\ntotal_reads,0\\n' > ${prefix}.genome_individual.csv
    """
}
