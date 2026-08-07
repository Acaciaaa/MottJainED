#!/usr/bin/env bash
# generator：建立/复用 point 的 ED 快照并拟合固定 Lambda。
#SBATCH --job-name=mj-generator
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=sdicnormal
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%x-%j.out

# ===== 用户配置区 =====
config="config/my_run.toml"
point="nm5_candidate_01"
threads=16
# =====================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
(( threads >= 1 && threads <= allocated_cpus )) || { echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2; exit 2; }
export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"
julia --startup-file=no --project=. bin/mottjain.jl generator \
  --config="$config" --point="$point"
