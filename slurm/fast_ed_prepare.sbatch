#!/usr/bin/env bash
#SBATCH --job-name=mj-n7-cache
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# N=7 first run requests 16 CPUs mainly for the partition's per-CPU memory.
# Computation remains at 8 threads until seff establishes a safe smaller request.
#SBATCH --cpus-per-task=16
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
julia --startup-file=no --project=. scripts/fast_ed.jl prepare --config="$CONFIG"
