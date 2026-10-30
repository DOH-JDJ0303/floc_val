# floc-consistency

Tests how consistent [floc](https://github.com/NW-PaGe/floc) clustering is when the
same assemblies are added in batches of different sizes, and measures SKA2 core SNP
distances within floc's final clusters and between samples floc groups inconsistently.

## How it works

1. The assemblies are put in one shuffled order (seeded), used for every batch size.
2. `FLOC_BATCHES` (one task per batch size, floc container) runs every floc iteration
   inside the task: batch 0 is clustered, each later batch is added to the previous
   batch's signatures (`floc -o out_i new_i/ out_{i-1}/sigs`). Only the final
   `clusters.csv` is kept. The run with all samples in one batch is the reference.
3. `PLAN_PAIRS` finds the pairs to measure: discrepant pairs (grouped differently from
   the reference in at least one run) and every pair within each run's final clusters.
4. `SKA_BUILD` (one task per sample, StaPH-B ska2) and `SKA_DISTANCE` (one task per
   core group) compute core SNP distances.
5. `SUMMARISE` writes the tables.

## Run

```bash
# locally with Docker
nextflow run main.nf --assemblies assemblies/ --outdir results --batch_sizes 10,50

# extra floc options
nextflow run main.nf --assemblies assemblies/ --outdir results --batch_sizes 10,50 \
    --floc_args '--min-hash-frac 0'

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
| `--batch_sizes` | | Comma-separated batch sizes, e.g. `10,50` (overrides start/step/stop) |
| `--start`, `--step`, `--stop` | 1, 10, all | Batch-size schedule when `--batch_sizes` is not set |
| `--shuffle`, `--seed` | true, 1 | Sample order |
| `--floc_args` | | Extra floc options for every call |
| `--core` | `group` | SKA2 core: `group`, `global` or `pairwise` |
| `--ska_k` | 31 | SKA2 split k-mer size |
| `--snp_threshold` | 100 | Count discrepant pairs with fewer core SNPs than this |
| `--floc_container` | `public.ecr.aws/o8h2f0o1/floc:1.2.1` | |
| `--ska_container` | `staphb/ska2:0.5.1` | |

## Outputs

| File | Contents |
|---|---|
| `summary.tsv` | n, min, max, mean, sd, median of core SNPs: discrepant pairs (all, together/apart in the reference) and within-cluster pairs (per batch size and over all runs). Header lines give counts, including discrepant pairs under `--snp_threshold`. |
| `discrepant_pairs.tsv` | One row per discrepant pair, with the batch sizes that differ from the reference |
| `within_cluster_pairs.tsv` | One row per pair in the same final cluster, per batch size |
| `cluster_stats.tsv` | Per batch size and cluster: size, core SNP statistics, members |
| `floc/bs<N>.clusters.csv`, `floc/bs<N>.floc.log` | Final floc clusters and floc log for each batch size |
| `sample_order.txt` | The sample order used |
