#!/usr/bin/env bash
# 本文件只用于采用 Slurm 的服务器。Mac 本地运行不需要它。

# 以下 #SBATCH 行是提交给 Slurm 的资源申请，不是 Julia 代码。
#SBATCH --job-name=mottjain-fss
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=mottjain-%j.out

# 任一命令失败就停止；使用未定义变量也停止，避免失败后继续写出假结果。
set -euo pipefail

# 无论从哪个目录 sbatch，都找到本脚本上一层的 MottJainED 根目录。
project_dir="$(cd "$(dirname "$0")/.." && pwd)"

# Julia/FuzzifiED 使用 Slurm 分配的 CPU；BLAS 保持单线程，避免线程过量。
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1

# 进入项目根目录，然后执行 fss-all。
# 正式使用时通常应把 default.toml 改成你自己的配置文件。
cd "$project_dir"
julia --startup-file=no --project=. bin/mottjain.jl fss-all \
  --config=config/default.toml --override=config/fss_profiles/fss7.toml
