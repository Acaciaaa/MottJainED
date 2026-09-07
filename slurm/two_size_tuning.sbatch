#!/usr/bin/env bash
# 独立的 N=5/6 matching 调参任务；N=3/4 只给 muc continuation 导航。
#SBATCH --job-name=mj-2size-uf0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-%x-%j.out

# ===== 用户配置区 =====
config="config/my_run.toml"
profile="config/two_size_tuning/uf0_after_vf0_stage2.toml"
threads=16
# =====================

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
