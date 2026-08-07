#!/usr/bin/env bash
# optimize：按指定 profile 优化 Hamiltonian 参数。
#SBATCH --job-name=mj-optimize
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=sdicnormal
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%x-%j.out

# ===== 用户配置区 =====
config="config/my_run.toml"
profile="config/optimization_profiles/mu_uf_uf0_vf0.toml"
nm1=5
k=10
# =====================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"
julia --startup-file=no --project=. bin/mottjain.jl optimize \
  --config="$config" --override="$profile" --nm1="$nm1" --k="$k"
