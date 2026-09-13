#!/usr/bin/env bash
# Run package loading, bounded cache retirement, and pipeline initialization on
# one allocated CPU instead of blocking the login node.
#SBATCH --job-name=mj-n7-bootstrap
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sdicnormal
#SBATCH --time=00:30:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
BASE_CONFIG="${BASE_CONFIG:?BASE_CONFIG is required}"
cd "$PROJECT_ROOT"

# Match precompilation to this job's allocation.
export JULIA_NUM_THREADS=1
export JULIA_NUM_PRECOMPILE_TASKS=1

julia --startup-file=no --project=. scripts/fast_ed_pipeline.jl \
    retire-configured-caches --config="$BASE_CONFIG"

bash scripts/submit_fast_ed_pipeline.sh "$BASE_CONFIG"
