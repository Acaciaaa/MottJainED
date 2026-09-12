#!/usr/bin/env bash
# One Hamiltonian only. Four processes each hold one sector's persistent matrices.
# 4 x 8 threads = 32 CPUs; the controller is idle during eigensolver calls.
#SBATCH --job-name=mj-n7-uf0165-search
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
CONFIG="config/fast_ed/n7_uf0_165_search.toml"
threads=8
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
(( allocated_cpus >= 4 * threads )) || {
    echo "This four-sector resident search requires 32 allocated CPUs" >&2
    exit 2
}
export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$PROJECT_ROOT"
julia --startup-file=no --project=. scripts/fast_mu_search.jl run \
    --config="$CONFIG" --threads-per-worker="$threads"
