#!/usr/bin/env bash
# Recover a reviewed pipeline only when its fitted vertex was snapped onto an
# already sampled point.  The normal controller then submits one fresh mu and
# its dependent final controller.
#SBATCH --job-name=mj-n7-recover
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sdicnormal
#SBATCH --time=00:30:00
#SBATCH --output=slurm-logs/%x-%j.out

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
BASE_CONFIG="${BASE_CONFIG:?BASE_CONFIG is required}"
cd "$PROJECT_ROOT"

julia --startup-file=no --project=. scripts/fast_ed_pipeline.jl \
    resume-stalled-fit --config="$BASE_CONFIG" >/dev/null

exec bash slurm/fast_ed_pipeline_controller.sbatch
