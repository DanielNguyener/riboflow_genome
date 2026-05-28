#  — riboflow_genome (DSL2)

Ribosome-profiling **genome** alignment pipeline, migrated from the DSL1
monolith `RiboFlow.groovy` to **Nextflow DSL2** in an nf-core-style layout.
`RiboFlow.groovy` is retained at the repo root **as a parity reference only** —
it is DSL1 and will NOT run under this env (Nextflow ≥24 / Java 17).

## Layout
- `main.nf` → `workflows/riboflow.nf` → `subworkflows/local/*` → `modules/local/*`
- `conf/base.config` (defaults), `conf/local.config` (executor budget + per-tool
  resources), `conf/modules.config` (per-module storeDir/publishDir/ext), `conf/test.config`
- `lib/Utils.groovy` — Groovy helpers ported from `RiboFlow.groovy:7-98`
- Docs: `` (DSL1→DSL2 map + status), ``
  (meta map + channel shapes), `docs/ARCHITECTURE.md`

## Conventions
- **One tool per module.** Each `modules/local/*.nf` has `tag`, typed I/O, a
  `script:` and a `stub:` block.
- **Canonical meta map:** per-lane `meta = [id, lane, strand:'F']`; per-sample
  `meta = [id, strand:'F']`. See ``.
- **CLI args / output names / dir locations are NOT hardcoded in modules** —
  they come from `conf/modules.config` via `ext.args` / `ext.prefix` and the
  storeDir/publishDir closures. Shared modules (genome + transcriptome) are
  parametrized by `ext` flags: `presort`, `append_index`, `emit_counts`,
  `emit_bed_counts`, `use_merged_counts`.
- **`lib/Utils` is NOT available inside config directive closures** (they run in
  the task context). `conf/modules.config` therefore inlines the storeDir/
  publishDir path logic using `params` only. `Utils.*` is fine in `.nf` scripts.
- **Parity rule:** keep all output/intermediate directory names and file names
  identical to the DSL1 pipeline (incl. `transcriptome_alignment`).

## Environment (single consolidated env)
- `environment.yaml` = ONE conda env: Nextflow ≥24, `openjdk>=17`, **umicollapse
  via bioconda** (no more hand-shipped jar / `java11` wrapper), all bio tools at
  latest-compatible versions. `docker/Dockerfile` builds purely from it.
- Profiles: `standard` (ambient env), `conda` (Nextflow manages env), `docker`,
  `test` (stub fixtures).

## Run
- Stub wiring check (no tools needed):
  `nextflow run main.nf -stub-run -profile test [--dedup_method none|position|umicollapse] [--star.output_transcriptome_bam true]`
- Real run: `nextflow run main.nf -profile conda -params-file example_position.yaml`

## Scope
- **Phase 1 (here): genome path only** + optional STAR transcriptome-*projected*
  BAM dedup (`do_tx_dedup`, BAM/BED only).
- **Deferred:** bowtie2 transcriptome alignment; `ribopy create` / `.ribo`;
  the RNA-seq path. RNA-seq params are accepted but ignored with a warning.
