#!/usr/bin/env bash
# Submit with --array=0-(4*number_of_mu_points-1)%4. The supplied five-mu scout
# therefore uses --array=0-19%4; the one-mu k=30 check uses --array=0-3%4.
#SBATCH --job-name=mj-n7-sector
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# Reserve extra CPUs for N=7 memory, but keep the measured eigensolver optimum
# at 8 threads. Revisit this after checking MaxRSS with seff.
#SBATCH --cpus-per-task=16
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-logs/%x-%A_%a.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
CONFIG="${CONFIG:-config/fast_ed/n7_retained_k20.toml}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?submit this script with sbatch --array=...}"
MU_INDEX=$((TASK_ID / 4 + 1))
SECTOR_INDEX=$((TASK_ID % 4 + 1))
threads="${THREADS:-8}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
(( threads >= 1 && threads <= allocated_cpus )) || {
    echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2
    exit 2
}
export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$PROJECT_ROOT"
julia --startup-file=no --project=. scripts/fast_ed.jl solve \
    --config="$CONFIG" \
    --mu-index="$MU_INDEX" \
    --sector-index="$SECTOR_INDEX"
