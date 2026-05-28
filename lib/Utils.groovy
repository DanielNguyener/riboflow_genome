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
}
