#!/usr/bin/env bash
#SBATCH --job-name=so3-n6-conf-check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-logs/so3-n6-conf-check-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
result_parent="output/so3lver/conformal_optimization"
result_name="n6_full_algebra_pilot_job-${SLURM_JOB_ID}"
result_dir="$result_parent/$result_name"
mkdir -p "$result_dir"

# Fixed anchor, old corrected-search winner, and a cold repeat of that winner.
# This does not start a Hamiltonian search and does not submit another job.
pilot_status=0
srun --ntasks=1 --cpus-per-task=1 julia --startup-file=no --project=. \
  scripts/so3lver_conformal_optimize.jl \
  --config=config/so3lver/n6_projected_conformal_optimization.toml \
  --output="$result_dir" --pilot-only=true \
  --pilot-point=0.5338034014111593,2.2231485242157642,0.24493901605947094,0.3223055538287333 \
  > "$result_dir/pilot.stdout.log" 2>&1 || pilot_status=$?
archive="$result_parent/${result_name}.tar.gz"
tar -czf "$archive" -C "$result_parent" "$result_name"
sha256sum "$archive" > "$archive.sha256"
echo "pilot_exit_status=$pilot_status"
echo "result_archive=$PWD/$archive"
echo "result_checksum=$PWD/$archive.sha256"
exit "$pilot_status"
