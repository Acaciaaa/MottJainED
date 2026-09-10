#!/usr/bin/env bash
# Observable check for the retained stage-8 point. The three independent
# tasks run scaling at N=6, gap at N=3,4,5,6, and density at N=6.
# Existing stage-8 data already supply FSS, so it is not recomputed here.
#SBATCH --job-name=mj-pointcheck
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%x-%A_%a.out

# ===== User settings =====
config="config/my_run.toml"
profile="config/evaluation_profiles/recommended_uf0_1834_vf0_055_v0_034.toml"
commands=("scaling" "gap" "density")
threads=16
# =========================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
task_index="${SLURM_ARRAY_TASK_ID:-0}"

(( task_index >= 0 && task_index < ${#commands[@]} )) || {
  echo "invalid array index: $task_index" >&2
  exit 2
}
command="${commands[$task_index]}"
(( threads >= 1 && threads <= allocated_cpus )) || {
  echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2
  exit 2
}

export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$project_dir"
julia --startup-file=no --project=. bin/mottjain.jl "$command" \
  --config="$config" --override="$profile"
