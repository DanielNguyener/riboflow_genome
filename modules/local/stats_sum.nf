// Sum per-lane stats into a per-sample row. Ports
// `sum_individual_genome_alignment_stats` (RiboFlow.groovy:1677-1710). When
// ext.use_merged_counts=true (dedup position/umicollapse) the merged dedup
// counts override the summed per-lane dedup rows via `rfc update-dedup-counts`;
// otherwise (dedup none) the placeholder count files are ignored.
process STATS_SUM {
    executor 'local'
    tag "${meta.id}"

    input:
    tuple val(meta), path(stat_files),
          path(dedup_total,     stageAs: 'merged_dedup.total.count'),
          path(dedup_primary,   stageAs: 'merged_dedup.primary.count'),
          path(dedup_secondary, stageAs: 'merged_dedup.secondary.count'),
          path(dedup_unique,    stageAs: 'merged_dedup.unique.count')

    output:
    tuple val(meta), path("${meta.id}.${task.ext.label ?: 'genome'}_merged.csv"), emit: csv

    script:
    def label       = task.ext.label ?: 'genome'
    def use_counts  = task.ext.use_merged_counts ?: false
    def unique_only = (task.ext.unique_only != null) ? task.ext.unique_only : true
    // In multi-mapper mode, patch dedup_unique_alignments and dedup_multi_primary_alignments
    // inline via Python so this works with any rfc version (old Docker images don't have
    // --dedup-unique-file in rfc update-dedup-counts).
    def update_unique_cmd = unique_only ? '' : """
python3 - << 'PYEOF'
import csv
fname = '${meta.id}.${label}_merged.csv'
unique = int(open('merged_dedup.unique.count').read().strip().split()[0])
rows = list(csv.reader(open(fname)))
primary = next((int(r[1]) for r in rows[1:] if r and r[0] == 'dedup_primary_alignments'), None)
for r in rows[1:]:
    if r and len(r) >= 2:
        if r[0] == 'dedup_unique_alignments':
            r[1] = str(unique)
        elif r[0] == 'dedup_multi_primary_alignments' and primary is not None:
            r[1] = str(primary - unique)
with open(fname, 'w', newline='') as f:
    csv.writer(f).writerows(rows)
PYEOF
"""
    if (use_counts) {
        """
        rfc sum-stats -n ${meta.id} -o ${meta.id}.${label}_merged.tmp.csv ${stat_files}
        rfc update-dedup-counts \\
            --dedup-total-file merged_dedup.total.count \\
            --dedup-primary-file merged_dedup.primary.count \\
            --dedup-secondary-file merged_dedup.secondary.count \\
            --input-csv ${meta.id}.${label}_merged.tmp.csv \\
            --output-csv ${meta.id}.${label}_merged.csv
        ${update_unique_cmd}
        """
    } else {
        """
        rfc sum-stats -n ${meta.id} -o ${meta.id}.${label}_merged.csv ${stat_files}
        """
    }

    stub:
    def label = task.ext.label ?: 'genome'
    """
    printf ',${meta.id}\\ntotal_reads,0\\n' > ${meta.id}.${label}_merged.csv
    """
}
