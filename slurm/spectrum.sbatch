#!/usr/bin/env bash
# spectrum：单个系统大小上的完整低能谱扫描。
#SBATCH --job-name=mj-spectrum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=sdicnormal
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out

# ===== 用户配置区：通常只改这里和上面的资源 =====
config="config/my_run.toml"
# =================================================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"
julia --startup-file=no --project=. bin/mottjain.jl spectrum --config="$config"
