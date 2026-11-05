#!/usr/bin/env nextflow
/*
 * floc-consistency
 *
 * Tests how consistent floc clustering is when the same assemblies are added
 * in batches of different sizes, and compares floc's clusters with pairwise
 * SNP distances.
 *
 *   - SNP distances between every pair of samples are computed first with
 *     SKA2 (SKA_BUILD + SKA_DISTANCE), or read from --distance_matrix
 *   - at most --max_samples assemblies (default 50) are picked at random
 *     (seeded) and used for every step
 *   - one shuffled sample order (seeded) is used for every batch size
 *   - FLOC_BATCHES runs all of floc's iterations for one batch size in a
 *     single task: batch 0 is clustered, each later batch is added to the
 *     previous batch's signatures; only the final clusters.csv is kept
 *   - the run with all samples in one batch is the reference
 *   - SUMMARISE reports intra- and inter-cluster SNP distances and pairs
 *     with < --snp_threshold SNPs that floc placed in different clusters
 *   - PLOT_TREE draws a neighbour-joining tree of the SNP distances with a
 *     column of cluster tiles per batch size (interactive plotly HTML)
 *
 * Usage:
 *   nextflow run main.nf --assemblies <dir> --outdir <dir> [--batch_sizes 10,50]
 *                        [--distance_matrix <file>]
 *                        [--max_samples 50]
 */
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
def sampleName(f) {
    // floc's naming: strip .gz, then one extension
    def b = f.name
    if (b.toLowerCase().endsWith('.gz')) b = b[0..-4]
    def i = b.lastIndexOf('.')
    return i > 0 ? b[0..<i] : b
}

def batchSizes(n) {
    def sizes = []
    if (params.batch_sizes) {
        sizes = params.batch_sizes.toString().tokenize(', ').collect { s -> Math.min(Math.max(1, s as int), n) }
    } else {
        def first = params.start as int
        def last = params.stop ? Math.min(params.stop as int, n) : n
        if (first <= last) sizes = (first..last).step(params.step as int)
    }
    sizes = sizes + [n]                 // all samples in one batch = reference
    return sizes.unique().sort()
}

// floc command-line options from the floc_* params. Value options that are unset
// (null) or set to 'none' are not passed, so floc uses its own default.
def flocArgs() {
    def opts = [
        floc_distance         : '--distance',
        floc_ksize            : '--ksize',
        floc_scaled           : '--scaled',
        floc_batch            : '--batch',
        floc_min_hash_freq    : '--min-hash-freq',
        floc_global_hash_freq : '--global-hash-freq',
        floc_min_hash_frac    : '--min-hash-frac',
        floc_assign_metric    : '--assign-metric',
        floc_assign_hash_freq : '--assign-hash-freq',
        floc_dbscan_metric    : '--dbscan-metric'
    ]
    def flags = [
        floc_ignore_qc        : '--ignore-qc',
        floc_skip_qc          : '--skip-qc',
        floc_overwrite        : '--overwrite',
        floc_plot             : '--plot',
        floc_pairwise         : '--pairwise'
    ]
    def args = []
    opts.each  { p, flag ->
        def v = params[p]
        if (v == null || v == '' || v.toString().toLowerCase() == 'none') return
        // Nextflow sets a parameter given without a value (or with '') to true
        if (v instanceof Boolean) error "--${p} needs a value (or 'none' to use floc's default)"
        args << "${flag} ${v}"
    }
    flags.each { p, flag -> if (params[p].toString().toLowerCase() == 'true') args << flag }
    return args.join(' ')
}

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------
process FLOC_BATCHES {
    tag "batch size ${batch_size}"
    label 'floc'
    publishDir "${params.outdir}/floc", mode: 'copy'

    input:
    val batch_size
    path order
    path assemblies, stageAs: 'assemblies/*'

    output:
    path "bs${batch_size}.clusters.csv", emit: clusters
    path "bs${batch_size}.floc.log",     emit: log

    script:
    """
    set -euo pipefail
    log=bs${batch_size}.floc.log
    : > "\$log"

    # sample name -> staged file (floc names samples by file name minus extensions)
    declare -A FILE
    for f in assemblies/*; do
        b=\$(basename "\$f"); b=\${b%.gz}; b=\${b%.*}
        FILE["\$b"]="\$PWD/\$f"
    done
    mapfile -t ORDER < ${order}
    n=\${#ORDER[@]}
    bs=${batch_size}

    i=0; prev=""
    for (( start=0; start<n; start+=bs )); do
        dir=\$(printf 'batch_%04d' \$i)
        mkdir -p "\$dir/new"
        for s in "\${ORDER[@]:start:bs}"; do
            ln -s "\${FILE[\$s]}" "\$dir/new/"
        done
        if [[ -z "\$prev" ]]; then
            inputs=( "\$dir/new" )
        else
            inputs=( "\$dir/new" "\$prev/out/sigs" )
        fi
        echo "### batch \$i: \$(ls "\$dir/new" | wc -l) new assemblies" >> "\$log"
        floc ${flocArgs()} -o "\$dir/out" "\${inputs[@]}" >> "\$log" 2>&1
        # floc exits 0 when a sample fails QC, so check its output as well
        if [[ ! -s "\$dir/out/clusters.csv" ]]; then
            echo "ERROR: floc wrote no clusters.csv for batch \$i (see \$log)" >&2
            tail -n 20 "\$log" >&2
            exit 1
        fi
        [[ -n "\$prev" ]] && rm -rf "\$prev"     # previous batch no longer needed
        prev=\$dir; i=\$((i + 1))
    done
    cp "\$prev/out/clusters.csv" bs${batch_size}.clusters.csv
    """
}

process SKA_BUILD {
    tag "${sample}"
    label 'ska'
    label 'small'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}.skf")

    script:
    """
    printf '%s\\t%s\\n' "${sample}" "${assembly}" > input.tsv
    ska build -f input.tsv -o "${sample}" -k ${params.ska_k} --threads ${task.cpus}
    """
}

process SKA_DISTANCE {
    tag "${skfs.size()} genomes"
    label 'ska'
    publishDir "${params.outdir}/distances", mode: 'copy'

    input:
    path skfs, stageAs: 'skf/*'

    output:
    path 'ska_distances.tsv'

    script:
    // global: split k-mers present in every sample; pairwise: those shared by each pair
    def min_freq = params.core == 'global' ? 1 : 0
    """
    ska merge -o all skf/*.skf
    ska distance --min-freq ${min_freq} --threads ${task.cpus} -o ska_distances.tsv all.skf
    rm -f all.skf
    """
}

process SUMMARISE {
    label 'floc'
    label 'small'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path clusters
    path distances
    val n_samples

    output:
    path 'summary.tsv'
    path 'run_metrics.tsv'
    path 'split_close_pairs.tsv'
    path 'discrepant_pairs.tsv'
    path 'within_cluster_pairs.tsv'
    path 'cluster_stats.tsv'

    script:
    def source = params.distance_matrix
        ? "supplied matrix ${file(params.distance_matrix).name}"
        : "SKA2 k=${params.ska_k}, core=${params.core}"
    """
    floc_consistency.py summarise --clusters ${clusters} --distances ${distances} \\
        --n-samples ${n_samples} --snp-threshold ${params.snp_threshold} \\
        --distance-source '${source}'
    """
}

process PLOT_TREE {
    label 'floc'
    label 'small'
    publishDir "${params.outdir}/tree", mode: 'copy'

    input:
    path clusters
    path distances
    val n_samples

    output:
    path 'nj_tree.nwk'
    path 'nj_tree_clusters.html'

    script:
    """
    floc_consistency.py tree --clusters ${clusters} --distances ${distances} \\
        --n-samples ${n_samples} --snp-threshold ${params.snp_threshold} \\
        --newick nj_tree.nwk --html nj_tree_clusters.html
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    if (!params.assemblies) error "Please provide --assemblies <directory of assemblies>"
    if (!(params.core in ['global', 'pairwise'])) error "--core must be global or pairwise"
    if (params.floc_assign_metric && !(params.floc_assign_metric in ['avg', 'min', 'max', 'query']))
        error "--floc_assign_metric must be avg, min, max or query"
    if (params.floc_dbscan_metric && !(params.floc_dbscan_metric in ['avg', 'min', 'max']))
        error "--floc_dbscan_metric must be avg, min or max"
    if (params.containsKey('floc_args')) error "--floc_args has been replaced by individual --floc_* options (see README)"
    log.info "floc-consistency: floc options: ${flocArgs() ?: '(floc defaults)'}"

    def dir = params.assemblies.toString().replaceAll('/+$', '')
    ch_asm = Channel.fromPath("${dir}/*.{fa,fasta,fna,fa.gz,fasta.gz,fna.gz}", checkIfExists: true)

    // random subset of at most --max_samples assemblies (seeded); 0 = all
    ch_all = ch_asm
        .toSortedList { a, b -> a.name <=> b.name }
        .map { files ->
            def names = files.collect { sampleName(it) }
            def dups = names.countBy { it }.findAll { it.value > 1 }.keySet()
            if (dups) error "Duplicate sample names (file name without extensions): ${dups.join(', ')}"
            def max = (params.max_samples ?: 0) as int
            if (max <= 0 || files.size() <= max) return files
            def picked = new ArrayList(files)
            Collections.shuffle(picked, new Random(params.seed as long))
            log.info "floc-consistency: randomly selected ${max} of ${files.size()} assemblies (--max_samples ${max}, --seed ${params.seed})"
            picked.take(max).sort { it.name }
        }

    // 1. pairwise SNP distances between every pair of samples
    if (params.distance_matrix) {
        log.info "floc-consistency: using the supplied distance matrix ${params.distance_matrix}; SKA2 is skipped"
        ch_dist = Channel.fromPath(params.distance_matrix, checkIfExists: true).first()
    } else {
        SKA_BUILD(ch_all.flatMap { it }.map { f -> [sampleName(f), f] })
        SKA_DISTANCE(SKA_BUILD.out.map { s, skf -> skf }.collect())
        ch_dist = SKA_DISTANCE.out
    }

    // 2. floc clustering at each batch size
    ch_design = ch_all.map { files ->
        def names = files.collect { sampleName(it) }
        def order = new ArrayList(names)
        if (params.shuffle) Collections.shuffle(order, new Random(params.seed as long))
        def n = files.size()
        log.info "floc-consistency: ${n} assemblies; batch sizes ${batchSizes(n).join(', ')}"
        [order, batchSizes(n), n]
    }

    ch_order = ch_design
        .map { it[0].join('\n') + '\n' }
        .collectFile(name: 'sample_order.txt', storeDir: "${params.outdir}")
        .first()
    ch_sizes = ch_design.flatMap { it[1] }
    ch_n     = ch_design.map { it[2] }

    // one task per batch size; floc iterations run inside the task
    FLOC_BATCHES(ch_sizes, ch_order, ch_all)
    ch_clusters = FLOC_BATCHES.out.clusters.collect()

    // 3. intra-/inter-cluster distances and close pairs split across clusters
    SUMMARISE(ch_clusters, ch_dist, ch_n)

    // 4. neighbour-joining tree of the SNP distances with each batch size's clusters
    PLOT_TREE(ch_clusters, ch_dist, ch_n)
}
