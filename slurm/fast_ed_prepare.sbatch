#!/usr/bin/env bash
#SBATCH --job-name=mj-n7-cache
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# Recorded N=7 prepare MaxRSS fits the partition's 8-CPU memory allowance.
#SBATCH --cpus-per-task=8
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
CONFIG="${CONFIG:-config/fast_ed/n7_retained_k20.toml}"
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
extra_args=(--config="$CONFIG")
case "${FORCE_REBUILD:-false}" in
    true) extra_args+=(--force) ;;
    false) ;;
    *) echo "FORCE_REBUILD must be true or false" >&2; exit 2 ;;
esac
echo "$(date -Is) prepare: config=$CONFIG threads=$threads rebuild=${FORCE_REBUILD:-false} host=$(hostname)"
julia --startup-file=no --project=. scripts/fast_ed.jl prepare "${extra_args[@]}"
