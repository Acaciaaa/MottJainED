#!/usr/bin/env bash
# S-only iterative FSS stage 8: Uf scan with U0=9*Uf.
# The three scan points share each size's Hamiltonian cache in one job.
# N=6 audits the current Uf=0.46 point against the wide muc search first.
#SBATCH --job-name=mj-s-stage8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-%x-%j.out

# ===== User settings =====
config="config/my_run.toml"
profile="config/two_size_tuning/s_stage8_uf_linked.toml"
threads=16
# =========================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"

(( threads >= 1 && threads <= allocated_cpus )) || {
  echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2
  exit 2
}

export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$project_dir"
julia --startup-file=no --project=. scripts/two_size_tuning.jl \
  --config="$config" --profile="$profile"
