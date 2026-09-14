#!/usr/bin/env bash
# Four independent Hamiltonians. Every task stops after its N=3,4 guides and
# N=5,6 matching run; deciding whether to run N=7 is a separate review step.
#SBATCH --job-name=mj-orig-n56
#SBATCH --array=0-3%4
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
    "config/two_size_tuning/n56_original_uf0_165.toml"
    "config/two_size_tuning/n56_original_uf0_200.toml"
    "config/two_size_tuning/n56_original_vf0_030.toml"
    "config/two_size_tuning/n56_original_vf0_052.toml"
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
