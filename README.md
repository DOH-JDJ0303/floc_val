# floc-consistency

Tests how consistent [floc](https://github.com/NW-PaGe/floc) clustering is when the
same assemblies are added in batches of different sizes, and checks floc's clusters
against pairwise SNP distances between every pair of samples.

## How it works

0. **Sample selection.** At most `--max_samples` (default 50) assemblies are picked at
   random (seeded with `--seed`); `--max_samples 0` uses all of them. Samples outside
   the selection are ignored everywhere, including rows of a supplied distance matrix.
1. **SNP distances.** `SKA_BUILD` (one task per sample, StaPH-B ska2) and
   `SKA_DISTANCE` (one task, all samples merged) compute SNP distances between every
   pair of samples. This runs alongside clustering. With `--distance_matrix`, the
   supplied matrix is used and SKA2 is skipped.
   Only the assemblies picked in step 0 are used, here and in every later step.
2. The assemblies are put in one shuffled order (seeded), used for every batch size.
3. `FLOC_BATCHES` (one task per batch size, floc container) runs every floc iteration
   inside the task: batch 0 is clustered, each later batch is added to the previous
   batch's signatures (`floc -o out_i new_i/ out_{i-1}/sigs`). Only the final
   `clusters.csv` is kept. The run with all samples in one batch is the reference.
4. `SUMMARISE` reports, for each batch size, SNP distances within clusters
   (intra-cluster) and between clusters (inter-cluster), and every pair of samples
   with fewer than `--snp_threshold` (default 100) SNPs that floc placed in
   different clusters.

## Run

```bash
# locally with Docker
nextflow run main.nf --assemblies assemblies/ --outdir results --batch_sizes 10,50

# extra floc options
nextflow run main.nf --assemblies assemblies/ --outdir results --batch_sizes 10,50 \
    --floc_args '--min-hash-frac 0'

# reuse distances from an earlier run (or any matrix) and skip SKA2
nextflow run main.nf --assemblies assemblies/ --outdir results2 --batch_sizes 10,50 \
    --distance_matrix results/distances/ska_distances.tsv

# HPC
nextflow run main.nf -profile apptainer --assemblies assemblies/ --outdir results
```

On Seqera Platform, add the repository as a pipeline, choose a cloud compute
environment, and set `assemblies` and `outdir` to bucket paths
(e.g. `s3://bucket/assemblies`, `s3://bucket/floc-consistency`). Every process
runs in a public container image, so no extra setup is needed.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--assemblies` | (required) | Directory of `.fa`/`.fasta`/`.fna` assemblies, optionally `.gz` |
| `--outdir` | `results` | Output directory |
| `--distance_matrix` | | All-vs-all SNP distances to use instead of SKA2 (see below) |
| `--max_samples` | 50 | Randomly select at most this many assemblies (seeded); `0` = all |
| `--batch_sizes` | | Comma-separated batch sizes, e.g. `10,50` (overrides start/step/stop) |
| `--start`, `--step`, `--stop` | 1, 10, all | Batch-size schedule when `--batch_sizes` is not set |
| `--shuffle`, `--seed` | true, 1 | Sample order |
| `--floc_args` | | Extra floc options for every call |
| `--core` | `pairwise` | SKA2 sites: `pairwise` (split k-mers shared by each pair) or `global` (present in every sample) |
| `--ska_k` | 31 | SKA2 split k-mer size |
| `--snp_threshold` | 100 | Pairs with fewer SNPs than this are flagged when floc places them in different clusters |
| `--floc_container` | `public.ecr.aws/o8h2f0o1/floc:1.2.1` | |
| `--ska_container` | `staphb/ska2:0.5.1` | |

### Distance matrix format

`--distance_matrix` accepts either form, tab-, comma- or space-separated, optionally
gzipped. Sample names must match the assembly file names without extensions.

- **Long:** one pair per line: `sample_a`, `sample_b`, `distance`. A column named
  `Distance`, `snps` or `snp_distance` is used if present, otherwise the third column.
  `ska distance` output (including `distances/ska_distances.tsv` from an earlier run)
  works as is.
- **Square:** a header row of sample names (first cell ignored), then one row per
  sample starting with its name.

Pairs missing from the matrix are left out of the statistics and counted in
`run_metrics.tsv`.

## Outputs

| File | Contents |
|---|---|
| `summary.tsv` | n, min, max, mean, sd, median of SNP distances for intra-cluster and inter-cluster pairs per batch size, and for discrepant pairs. Header lines give, per batch size, the number of pairs with < `--snp_threshold` SNPs placed in different clusters. |
| `run_metrics.tsv` | One row per batch size: clusters, singletons, intra- and inter-cluster distance n/min/median/mean/max, pairs under the threshold, how many of them floc split (count and fraction), samples involved, pairs without a distance |
| `split_close_pairs.tsv` | Every pair with < `--snp_threshold` SNPs that floc placed in different clusters, per batch size, with both cluster IDs and whether the reference also split it |
| `cluster_stats.tsv` | Per batch size and cluster: size, intra-cluster SNP statistics, nearest sample in another cluster (SNPs, cluster, pair), members |
| `discrepant_pairs.tsv` | One row per pair grouped differently from the reference, with the batch sizes that differ and the pair's SNP distance |
| `within_cluster_pairs.tsv` | One row per pair in the same final cluster, per batch size |
| `distances/ska_distances.tsv` | SKA2 all-vs-all distances (not written when `--distance_matrix` is given); reusable as `--distance_matrix` |
| `floc/bs<N>.clusters.csv`, `floc/bs<N>.floc.log` | Final floc clusters and floc log for each batch size |
| `sample_order.txt` | The selected samples, in the order used |
