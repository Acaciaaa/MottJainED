#!/usr/bin/env bash
# spectrum：单个系统大小上的完整低能谱扫描。
#SBATCH --job-name=mj-spectrum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=sdicnormal
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out

# ===== 用户配置区：通常只改这里和上面的资源 =====
config="config/my_run.toml"
threads=16  # FuzzifiED 实际线程数；可以小于为取得内存而申请的 CPU 数
# =================================================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
(( threads >= 1 && threads <= allocated_cpus )) || { echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2; exit 2; }
export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"
julia --startup-file=no --project=. bin/mottjain.jl spectrum --config="$config"
