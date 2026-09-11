#!/usr/bin/env bash
#SBATCH --job-name=mj-n7-collect
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sdicnormal
#SBATCH --time=00:30:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
CONFIG="${CONFIG:-config/fast_ed/n7_retained_k20.toml}"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$PROJECT_ROOT"
julia --startup-file=no --project=. scripts/fast_ed.jl collect --config="$CONFIG"
