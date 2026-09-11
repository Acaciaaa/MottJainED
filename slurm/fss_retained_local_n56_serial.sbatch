#!/usr/bin/env bash
# First, low-storage stage of the retained-point local FSS campaign.
# Uf0 completes and passes its CSV audit before the Vf0 family begins.
#SBATCH --job-name=mj-fss-n56-serial
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=sdicnormal
#SBATCH --time=48:00:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
threads="${THREADS:-8}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
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

uf0_profile="config/fss_profiles/n56_retained_local_uf0.toml"
uf0_results="output/fast_ed/fss_n56_uf0/fss/fss_optimize_results.csv"
julia --startup-file=no --project=. bin/mottjain.jl fss-all \
    --config=config/my_run.toml --override="$uf0_profile" --method=optimize
julia --startup-file=no --project=. scripts/audit_fss_csv.jl \
    --config="$uf0_profile" --source="$uf0_results"

vf0_profile="config/fss_profiles/n56_retained_local_vf0.toml"
vf0_results="output/fast_ed/fss_n56_vf0/fss/fss_optimize_results.csv"
julia --startup-file=no --project=. bin/mottjain.jl fss-all \
    --config=config/my_run.toml --override="$vf0_profile" --method=optimize
julia --startup-file=no --project=. scripts/audit_fss_csv.jl \
    --config="$vf0_profile" --source="$vf0_results"
