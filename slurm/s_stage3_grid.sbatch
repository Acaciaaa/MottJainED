#!/usr/bin/env bash
# S-only iterative FSS stage 3: three Uf0 slices run as independent array tasks.
# Each slice scans V0; N=3,4 guide muc and only N=5,6 enter matching.
#SBATCH --job-name=mj-s-stage3
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-%x-%A_%a.out

# ===== User settings =====
config="config/my_run.toml"
profiles=(
  "config/two_size_tuning/s_stage3_uf0_150.toml"
  "config/two_size_tuning/s_stage3_uf0_1834.toml"
  "config/two_size_tuning/s_stage3_uf0_220.toml"
)
threads=16
# =========================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
task_index="${SLURM_ARRAY_TASK_ID:-0}"

(( task_index >= 0 && task_index < ${#profiles[@]} )) || {
  echo "invalid array index: $task_index" >&2
  exit 2
}
profile="${profiles[$task_index]}"
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
