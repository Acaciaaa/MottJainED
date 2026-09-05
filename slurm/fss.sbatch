#!/usr/bin/env bash
# FSS：用作业数组分别计算 Uf0/Vf0 两种一维扫描，随后各自画图和拟合。
#SBATCH --job-name=mj-fss
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-%x-%A_%a.out

# ===== 用户配置区 =====
config="config/my_run.toml"
profiles=("config/fss_profiles/no_w_uf0.toml" "config/fss_profiles/no_w_vf0.toml")
method="optimize"  # 改成 grid 时，每个外层点会硬算 my_run.toml 中的 41 个 mu
threads=16         # 56 CPU 取得内存，但 FuzzifiED 先只开 16 线程
# =====================

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
task_index="${SLURM_ARRAY_TASK_ID:-0}"
(( task_index >= 0 && task_index < ${#profiles[@]} )) || { echo "invalid array index: $task_index" >&2; exit 2; }
profile="${profiles[$task_index]}"
(( threads >= 1 && threads <= allocated_cpus )) || { echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2; exit 2; }
export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"
julia --startup-file=no --project=. bin/mottjain.jl fss-all \
  --config="$config" --override="$profile" --method="$method"
