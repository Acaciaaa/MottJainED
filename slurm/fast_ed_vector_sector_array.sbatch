#!/usr/bin/env bash
# Four independent standard sectors, each retaining k eigenvectors.
#SBATCH --job-name=mj-n7-k40vec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-logs/%x-%A_%a.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
CONFIG="${CONFIG:-config/fast_ed/n7_original_tower_k40.toml}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?submit this script with sbatch --array=0-3%4}"
SECTOR_INDEX=$((TASK_ID + 1))
threads="${THREADS:-8}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
(( TASK_ID >= 0 && TASK_ID < 4 )) || {
    echo "array task must be in 0:3" >&2
    exit 2
}
(( threads >= 1 && threads <= allocated_cpus )) || {
    echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2
    exit 2
}

export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$PROJECT_ROOT"
extra_args=(--config="$CONFIG" --sector-index="$SECTOR_INDEX")
case "${FORCE_SOLVE:-false}" in
    true) extra_args+=(--force) ;;
    false) ;;
    *) echo "FORCE_SOLVE must be true or false" >&2; exit 2 ;;
esac
echo "$(date -Is) vector solve: sector_index=$SECTOR_INDEX config=$CONFIG threads=$threads force=${FORCE_SOLVE:-false} host=$(hostname)"
julia --startup-file=no --project=. scripts/fast_ed.jl vector-solve "${extra_args[@]}"
