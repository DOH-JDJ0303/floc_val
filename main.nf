#!/usr/bin/env nextflow
/*
 * floc-consistency
 *
 * Tests how consistent floc clustering is when the same assemblies are added
 * in batches of different sizes, and measures SKA2 core SNP distances within
 * the final clusters and between samples floc groups inconsistently.
 *
 *   - one shuffled sample order (seeded) is used for every batch size
 *   - FLOC_BATCHES runs all of floc's iterations for one batch size in a
 *     single task: batch 0 is clustered, each later batch is added to the
 *     previous batch's signatures; only the final clusters.csv is kept
 *   - the run with all samples in one batch is the reference
 *
 * Usage:
 *   nextflow run main.nf --assemblies <dir> --outdir <dir> [--batch_sizes 10,50]
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
        floc ${params.floc_args} -o "\$dir/out" "\${inputs[@]}" >> "\$log" 2>&1
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

process PLAN_PAIRS {
    label 'floc'
    label 'small'

    input:
    path clusters
    val n_samples

    output:
    path 'ska_groups.tsv'

    script:
    """
    floc_consistency.py plan --clusters ${clusters} --n-samples ${n_samples} \\
        --core ${params.core} --out ska_groups.tsv
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
    tag "${group} (${skfs.size()} genomes)"
    label 'ska'

    input:
    tuple val(group), path(skfs, stageAs: 'skf/*')

    output:
    path "dist_${group}.tsv"

    script:
    def min_freq = params.core == 'pairwise' ? 0 : 1
    """
    ska merge -o merged skf/*.skf
    ska distance --min-freq ${min_freq} --threads ${task.cpus} -o dist_${group}.tsv merged.skf
    """
}

process SUMMARISE {
    label 'floc'
    label 'small'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path clusters
    path distances, stageAs: 'dist/*'
    val n_samples

    output:
    path 'summary.tsv'
    path 'discrepant_pairs.tsv'
    path 'within_cluster_pairs.tsv'
    path 'cluster_stats.tsv'

    script:
    def dists = distances ? "--distances dist/*" : ""
    """
    floc_consistency.py summarise --clusters ${clusters} ${dists} --n-samples ${n_samples} \\
        --snp-threshold ${params.snp_threshold} --core ${params.core} --ska-k ${params.ska_k}
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    if (!params.assemblies) error "Please provide --assemblies <directory of assemblies>"
    if (!(params.core in ['group', 'global', 'pairwise'])) error "--core must be group, global or pairwise"

    def dir = params.assemblies.toString().replaceAll('/+$', '')
    ch_asm = Channel.fromPath("${dir}/*.{fa,fasta,fna,fa.gz,fasta.gz,fna.gz}", checkIfExists: true)

    ch_all = ch_asm.toSortedList { a, b -> a.name <=> b.name }

    ch_design = ch_all.map { files ->
        def names = files.collect { sampleName(it) }
        def dups = names.countBy { it }.findAll { it.value > 1 }.keySet()
        if (dups) error "Duplicate sample names (file name without extensions): ${dups.join(', ')}"
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

    // pairs to measure, grouped for SKA
    PLAN_PAIRS(ch_clusters, ch_n)
    ch_members = PLAN_PAIRS.out
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.sample, row.group] }

    // split k-mer file for each sample that is measured
    ch_asm_by_name = ch_asm.map { f -> [sampleName(f), f] }
    ch_to_build = ch_members.map { s, g -> s }.unique().join(ch_asm_by_name)
    SKA_BUILD(ch_to_build)

    // core SNP distances per group
    ch_groups = ch_members
        .combine(SKA_BUILD.out, by: 0)            // [sample, group, skf]
        .map { s, g, skf -> [g, skf] }
        .groupTuple()
    SKA_DISTANCE(ch_groups)

    SUMMARISE(ch_clusters, SKA_DISTANCE.out.collect().ifEmpty([]), ch_n)
}
