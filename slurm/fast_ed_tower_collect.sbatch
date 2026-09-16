#!/usr/bin/env bash
# Assemble the four vector checkpoints, fit Lambda, and write tower overlaps.
# The snapshot is saved before fitting/analysis, so an analysis error does not
# require repeating the expensive diagonalizations.
#SBATCH --job-name=mj-n7-tower
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=sdicnormal
#SBATCH --time=24:00:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
CONFIG="${CONFIG:-config/fast_ed/n7_original_tower_k40.toml}"
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
echo "$(date -Is) assembling k=40 snapshot: config=$CONFIG threads=$threads host=$(hostname)"
julia --startup-file=no --project=. scripts/fast_ed.jl tower-assemble --config="$CONFIG"
echo "$(date -Is) running standard-family tower analysis"
julia --startup-file=no --project=. scripts/fast_ed.jl tower-analyze --config="$CONFIG"
