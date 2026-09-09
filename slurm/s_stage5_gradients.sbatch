#!/usr/bin/env bash
# S-only iterative FSS stage 5: measure Uf0 and V0 directions around two low-q cuts.
# Each profile is independent; N=3,4 guide muc and only N=5,6 enter matching.
#SBATCH --job-name=mj-s-stage5
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-%x-%A_%a.out

# ===== User settings =====
config="config/my_run.toml"
profiles=(
  "config/two_size_tuning/s_stage5_uf0_vf0_045_v0_110.toml"
  "config/two_size_tuning/s_stage5_uf0_vf0_060_v0_110.toml"
  "config/two_size_tuning/s_stage5_v0_uf0_1834_vf0_045.toml"
  "config/two_size_tuning/s_stage5_v0_uf0_1834_vf0_060.toml"
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
