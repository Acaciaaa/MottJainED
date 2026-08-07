#!/usr/bin/env bash
# FSS：计算 grid/optimize 数据，随后画图并尝试联合拟合。
#SBATCH --job-name=mj-fss
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-%x-%j.out

# ===== 用户配置区 =====
config="config/my_run.toml"
profile="config/fss_profiles/fss5.toml"  # 换七项时写 fss7.toml
method="both"                             # grid、optimize 或 both
threads=16                                 # 56 CPU 取得内存，但 FuzzifiED 先只开 16 线程
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
julia --startup-file=no --project=. bin/mottjain.jl fss-all \
  --config="$config" --override="$profile" --method="$method"
