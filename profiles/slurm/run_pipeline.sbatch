#!/bin/bash
# Runs Snakemake as a small, long SLURM job; Snakemake then submits every pipeline step as its own job
# (see profiles/slurm/config.yaml). Set the partition, account and environment below, then submit from
# the repository root:
#   sbatch profiles/slurm/run_pipeline.sbatch
#SBATCH --job-name=edqtl_controller
#SBATCH --output=logs/%x-%j.log
#SBATCH --partition=CPU_PARTITION
#SBATCH --account=ACCOUNT
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00

set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
eval "$(conda shell.bash hook)"
conda activate rna_editing_env      # any environment with Snakemake 9 + snakemake-executor-plugin-slurm

snakemake --profile profiles/slurm
