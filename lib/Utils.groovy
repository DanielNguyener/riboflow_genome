// Groovy helpers ported verbatim (behavior-preserving) from RiboFlow.groovy:7-98.
// These are static so they can be called from conf/modules.config closures and
// from the workflow scripts. The DSL1 versions read the `params` binding
// directly; in DSL2 we pass `params` in explicitly because library classes do
// not see the script binding.

class Utils {

    // get_storedir / get_publishdir — resolve an absolute intermediates/output
    // directory for a given logical output_type, honouring user overrides under
    // params.output.{intermediates,output}. (RiboFlow.groovy:7-37)
    static String storedir(Map params, String output_type, boolean is_rnaseq = false) {
        def out = (params.output ?: [:])
        def inter = (out.intermediates ?: [:])
        def base_path = inter.base ?: 'intermediates'
        if (is_rnaseq) {
            base_path = base_path + '/rnaseq'
        }
        return new File(base_path, (inter[output_type] ?: output_type)).getCanonicalPath()
    }

    static String publishdir(Map params, String output_type, boolean is_rnaseq = false) {
        def out = (params.output ?: [:])
        def outout = (out.output ?: [:])
        def base_path = outout.base ?: 'output'
        if (is_rnaseq) {
            base_path = base_path + '/rnaseq'
        }
        return new File(base_path, (outout[output_type] ?: output_type)).getCanonicalPath()
    }

    static String individual_dir(Map params) {
        return (params.output ?: [:]).individual_lane_directory ?: 'individual'
    }

    static String merged_dir(Map params) {
        return (params.output ?: [:]).merged_lane_directory ?: 'merged'
    }

    // get_dedup_method — normalise the dedup_method param, honouring the legacy
    // boolean `deduplicate` flag. (RiboFlow.groovy:39-62)
    static String dedup_method(String dedup_arg, String dedup_old) {
        def valid_methods = ['position', 'umicollapse', 'none']
        def dedup_param = dedup_arg.toLowerCase()

        if (dedup_param != 'none') {
            if (dedup_param in valid_methods) {
                return dedup_param
            } else {
                System.err.println('Invalid deduplication method ' + dedup_param +
                                   ' . Valid methods are: ' + valid_methods.join(','))
                System.exit(1)
            }
        } else {
            if (dedup_old.toLowerCase() != 'false') {
                return 'position'
            } else {
                return 'none'
            }
        }
    }

    static String resolve_dedup_method(Map params) {
        return dedup_method(
            (params.get('dedup_method', 'none')).toString(),
            (params.get('deduplicate', false)).toString())
    }

    static String resolve_rnaseq_dedup_method(Map params) {
        return dedup_method(
            (params.rnaseq?.get('dedup_method', 'none') ?: 'none').toString(),
            'false')
    }

    static boolean do_tx_dedup(Map params) {
        def m = resolve_dedup_method(params)
        return (m == 'umicollapse' || m == 'position') &&
               ((params.star ?: [:]).output_transcriptome_bam ?: false)
    }

    // Per-thread memory budget for `samtools sort`. (RiboFlow.groovy:64-79)
    static int samtools_sort_mem_per_thread_mb(task) {
        int sort_threads = Math.min(task.cpus as int, 8)
        int est = (int) (task.memory.toMega() * 0.7 / sort_threads)
        return Math.min(768, Math.max(64, est))
    }

    // Genome mapping-quality cutoff, null-safe (default 255). Use this instead of
    // `params.genome?.mapping_quality_cutoff as int ?: 255`: the Elvis operator
    // treats a cutoff of 0 (keep multimappers) as falsy and wrongly substitutes
    // 255, which flips `genome_unique_only` to true and breaks the dedup
    // unique/multi stats split.
    static int genome_mapq(Map params) {
        def v = (params.genome ?: [:]).mapping_quality_cutoff
        return (v != null) ? (v as int) : 255
    }

    // True when the genome path keeps only unique (MAPQ-255) alignments, so the
    // unique/multi/secondary stats breakdown is degenerate and can be synthesised.
    static boolean genome_unique_only(Map params) {
        return genome_mapq(params) >= 255
    }

    // RNA-seq genome mapping-quality cutoff, null-safe (default 4). Honours both the
    // nested (rnaseq.genome.mapping_quality_cutoff) and flat (rnaseq.mapping_quality_cutoff)
    // YAML shapes. Same Elvis-on-zero trap as genome_mapq: a deliberate cutoff of 0
    // (keep multimappers) must NOT fall through to the default.
    static int rnaseq_genome_mapq(Map params) {
        def rg = (params.rnaseq ?: [:])
        def v = (rg.genome ?: [:]).mapping_quality_cutoff
        if (v == null) v = rg.mapping_quality_cutoff
        return (v != null) ? (v as int) : 4
    }

    // RNA-seq genome unique-only mode: any positive cutoff keeps only unique reads
    // (cutoff 0 = keep multimappers → emit the unique/multi/secondary breakdown).
    static boolean rnaseq_genome_unique_only(Map params) {
        return rnaseq_genome_mapq(params) > 0
    }

    // ════════════════════════════════════════════════════════════════════════
    //  samtools_filter_arguments — post-alignment `samtools view` record filter
    //
    //  One string per alignment route, replacing the old `ribo_filter_flags` /
    //  `filter_flags` integers (which could only ever express a single -F mask).
    //
    //  NOTE: conf/modules.config carries a mirrored copy of these resolvers as
    //  file-local closures. That duplication is forced: config directive closures
    //  cannot see lib/ — a bare `Utils` there resolves to a groovy.util.ConfigObject.
    //  Keep the two in sync; process script blocks CAN call this class.
    // ════════════════════════════════════════════════════════════════════════

    // Default mask: unmapped(4) + secondary(256) + supplementary(2048).
    static final String DEFAULT_SAMTOOLS_FILTER_ARGS = '-F 2308'

    // `samtools view` options that select which RECORDS are kept, split by whether
    // they consume a following value. Anything outside these lists is rejected:
    // the pipeline owns the output stream (-b/-h/-o/-@/…), and -q is owned by
    // `<route>.mapping_quality_cutoff` because that same integer ALSO derives the
    // stats `unique_only` mode on the two genome routes (genome_unique_only /
    // rnaseq_genome_unique_only above). A -q here would let the actual filter and
    // the stats CSV shape disagree with no way to tell.
    static final List<String> SAMTOOLS_FILTER_OPTS_WITH_VALUE = [
        '-f', '--require-flags',
        '-F', '--excl-flags', '--exclude-flags',
        '-G',
        '-e', '--expr',
        '-d', '--tag',
        '-D', '--tag-file',
        '-L', '--target-file', '--region-file',
        '-r', '--read-group',
        '-R', '--read-group-file',
        '-N', '--qname-file',
        '-s', '--subsample', '--subsample-seed',
        '-m', '--min-qlen',
    ]

    static final List<String> SAMTOOLS_FILTER_FLAGS_NO_VALUE = [
        '-M', '--use-index',
        '-P', '--fetch-pairs',
    ]

    static final List<String> SAMTOOLS_MAPQ_OPTS = ['-q', '--min-MQ', '--min-mq']

    // Split a user argument string into tokens, honouring single and double quotes,
    // so a quoted filter expression (-e 'qlen>=25') survives as ONE token.
    static List<String> shell_split(String s) {
        def toks = []
        def cur = new StringBuilder()
        boolean started = false
        int i = 0
        int n = s.length()
        while (i < n) {
            String c = s.substring(i, i + 1)
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                if (started) { toks << cur.toString(); cur = new StringBuilder(); started = false }
                i++
            }
            else if (c == "'" || c == '"') {
                int j = s.indexOf(c, i + 1)
                if (j < 0) throw new IllegalArgumentException("unterminated ${c} quote")
                cur.append(s.substring(i + 1, j)); started = true; i = j + 1
            }
            else { cur.append(c); started = true; i++ }
        }
        if (started) toks << cur.toString()
        return toks
    }

    // Re-quote every token for safe interpolation into a process script. Called from
    // the module (NOT the config), so nothing a user writes can reach the shell as
    // syntax: `-F 2308; echo x` becomes '-F' '2308;' 'echo' 'x', which samtools
    // rejects as a bad argument rather than bash running it.
    static String shell_quote_args(String s) {
        if (s == null || s.trim().isEmpty()) return ''
        return shell_split(s).collect { "'" + it.replace("'", "'\\''") + "'" }.join(' ')
    }

    // Returns a list of human-readable problems; empty means valid. The caller
    // (workflows/riboflow.nf) turns a non-empty list into a startup `error`, so a bad
    // string fails in seconds instead of inside a task.
    static List<String> validate_samtools_filter_arguments(String args, String route) {
        def errors = []
        if (args == null) return errors
        def key = "${route}.samtools_filter_arguments"
        // `--route.samtools_filter_arguments ''` on the command line: Nextflow parses a
        // bare empty value as a boolean flag, so the string arrives as "true". Caught
        // by the stray-argument rule below anyway, but the cause is worth naming.
        if (args.toString() == 'true' || args.toString() == 'false') {
            return ["${key}: resolved to the boolean `${args}`. A bare `--${key} ''` on the " +
                    'command line is read by Nextflow as a flag, not an empty string — set the key ' +
                    'in your params YAML instead if you meant "no extra filtering".'].collect { it.toString() }
        }
        List<String> toks
        try {
            toks = shell_split(args.toString())
        }
        catch (IllegalArgumentException e) {
            return ["${key}: ${e.message} in ${args.inspect()}".toString()]
        }
        boolean expect_value = false
        toks.each { t ->
            if (expect_value) { expect_value = false; return }
            if (!t.startsWith('-')) {
                // samtools view's positional arguments are the input file and region
                // specs — both supplied by the pipeline. A stray one here means the
                // string is malformed, and it is how an injection attempt shows up
                // (`-F 2308; echo PWNED` leaves `echo` and `PWNED` dangling).
                errors << ("${key}: unexpected argument `${t}` — only options are allowed; RiboFlow " +
                           'supplies the input BAM and any region arguments itself.').toString()
                return
            }
            String name = t.contains('=') ? t.substring(0, t.indexOf('=')) : t
            if (SAMTOOLS_MAPQ_OPTS.contains(name)) {
                errors << ("${key}: `${name}` is not allowed here — set `${route}.mapping_quality_cutoff` " +
                           'instead. It is also what derives the unique-only stats mode, so the two must ' +
                           'not be able to disagree.').toString()
            }
            else if (SAMTOOLS_FILTER_OPTS_WITH_VALUE.contains(name)) {
                if (!t.contains('=')) expect_value = true
            }
            else if (!SAMTOOLS_FILTER_FLAGS_NO_VALUE.contains(name)) {
                errors << ("${key}: `${name}` is not an accepted read-filtering option. RiboFlow owns the " +
                           'output stream (-b/-h/-o/-@/-c/…); only record-selection options are allowed: ' +
                           (SAMTOOLS_FILTER_OPTS_WITH_VALUE + SAMTOOLS_FILTER_FLAGS_NO_VALUE).join(' ')).toString()
            }
        }
        if (expect_value) {
            errors << "${key}: `${toks[-1]}` expects a value but none follows.".toString()
        }
        return errors
    }

    // The removed integer keys. Fail loudly with the translation rather than letting a
    // stale YAML silently fall through to the default mask.
    static List<String> legacy_filter_flag_errors(Map params) {
        def errors = []
        def check = { Map holder, String old_key, String old_path, String new_path ->
            if (holder != null && holder.containsKey(old_key)) {
                errors << ("`${old_path}` was removed. Use `${new_path}: \"-F ${holder[old_key]}\"` " +
                           'instead — it takes any samtools view read-filtering arguments, not just a ' +
                           'single -F mask.').toString()
            }
        }
        def rnaseq = (params.rnaseq ?: [:]) as Map
        check((params.genome ?: [:]) as Map, 'ribo_filter_flags',
              'genome.ribo_filter_flags', 'genome.samtools_filter_arguments')
        check((rnaseq.genome ?: [:]) as Map, 'filter_flags',
              'rnaseq.genome.filter_flags', 'rnaseq.genome.samtools_filter_arguments')
        check(rnaseq, 'filter_flags',
              'rnaseq.filter_flags', 'rnaseq.genome.samtools_filter_arguments')
        return errors
    }

    // Per-route resolvers. `!= null`, NEVER Elvis: an empty String is falsy in Groovy,
    // so `samtools_filter_arguments: ""` (deliberately filter nothing) must not fall
    // through to the default. Same trap as the Elvis-on-zero bug in genome_mapq above.
    static String genome_filter_args(Map params) {
        def v = (params.genome ?: [:]).samtools_filter_arguments
        return (v != null) ? v.toString() : DEFAULT_SAMTOOLS_FILTER_ARGS
    }

    static String transcriptome_filter_args(Map params) {
        def v = (params.transcriptome ?: [:]).samtools_filter_arguments
        return (v != null) ? v.toString() : DEFAULT_SAMTOOLS_FILTER_ARGS
    }

    // Honours both the nested (rnaseq.genome.*) and flat (rnaseq.*) YAML shapes,
    // matching rnaseq_genome_mapq.
    static String rnaseq_genome_filter_args(Map params) {
        def rg = (params.rnaseq ?: [:])
        def v = (rg.genome ?: [:]).samtools_filter_arguments
        if (v == null) v = rg.samtools_filter_arguments
        return (v != null) ? v.toString() : DEFAULT_SAMTOOLS_FILTER_ARGS
    }

    static String rnaseq_transcriptome_filter_args(Map params) {
        def v = ((params.rnaseq ?: [:]).transcriptome ?: [:]).samtools_filter_arguments
        return (v != null) ? v.toString() : DEFAULT_SAMTOOLS_FILTER_ARGS
    }
}
