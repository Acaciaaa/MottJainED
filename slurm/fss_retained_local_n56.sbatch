#!/usr/bin/env bash
# Reuse the stage-12 two-size workflow, including N=3,4 muc guides.
# Independent Uf0/Vf0 tasks; only subsequent N=7 points must run serially.
#SBATCH --job-name=mj-fss-n56-guided
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-logs/%x-%A_%a.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
threads="${THREADS:-8}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
profiles=(
    "config/two_size_tuning/n56_retained_local_uf0.toml"
    "config/two_size_tuning/n56_retained_local_vf0.toml"
)
task_index="${SLURM_ARRAY_TASK_ID:-0}"
(( task_index >= 0 && task_index < ${#profiles[@]} )) || {
    echo "invalid array index: $task_index" >&2
    exit 2
}
(( threads >= 1 && threads <= allocated_cpus )) || {
    echo "threads must be between 1 and allocated CPUs ($allocated_cpus)" >&2
    exit 2
}
export JULIA_NUM_THREADS="$threads"
export OMP_NUM_THREADS="$threads"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$PROJECT_ROOT"
mkdir -p slurm-logs

julia --startup-file=no --project=. scripts/two_size_tuning.jl \
    --config=config/my_run.toml --profile="${profiles[$task_index]}" --audit=true
