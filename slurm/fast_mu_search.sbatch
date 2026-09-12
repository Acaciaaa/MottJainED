#!/usr/bin/env bash
# One Hamiltonian and one eigensolver at a time; persistent matrices stay on disk.
# Reuse the measured 8-thread sector setting without reserving four idle workers.
#SBATCH --job-name=mj-n7-uf0165-search
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
CONFIG="config/fast_ed/n7_uf0_165_search.toml"
threads=8
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
(( allocated_cpus == threads )) || {
    echo "This sequential search requests exactly 8 CPUs and uses 8 threads" >&2
    exit 2
}
export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$PROJECT_ROOT"
julia --startup-file=no --project=. scripts/fast_mu_search.jl run \
    --config="$CONFIG"
