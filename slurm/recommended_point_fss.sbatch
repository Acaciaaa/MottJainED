#!/usr/bin/env bash
# Full FSS around the retained point: Uf0 and Vf0 scans run as two parallel jobs.
# Every scan value computes and plots N=3,4,5,6; no finite-size fit is performed.
#SBATCH --job-name=mj-rec-fss
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-%x-%A_%a.out

# ===== User settings =====
config="config/my_run.toml"
profiles=(
  "config/fss_profiles/recommended_point_uf0.toml"
  "config/fss_profiles/recommended_point_vf0.toml"
)
method="optimize"
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
julia --startup-file=no --project=. bin/mottjain.jl fss-all \
  --config="$config" --override="$profile" --method="$method"
