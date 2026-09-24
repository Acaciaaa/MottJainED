#!/usr/bin/env bash
#SBATCH --job-name=so3-n6-conf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-logs/so3-n6-conf-%j.out

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
export JULIA_NUM_THREADS="$allocated_cpus"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"

config_path="config/so3lver/n6_projected_conformal_optimization.toml"
result_parent="output/so3lver/conformal_optimization"
result_name="n6_multistart_01"
result_dir="$result_parent/$result_name"
mkdir -p "$result_dir"

# The driver resumes all completed Hamiltonian points from its signed CSV
# checkpoint.  Resubmitting this same file after timeout or preemption is safe.
julia --startup-file=no --project=. \
  scripts/so3lver_conformal_optimize.jl \
  --config="$config_path"

archive="$result_parent/${result_name}_job-${SLURM_JOB_ID}.tar.gz"
tar -czf "$archive" -C "$result_parent" "$result_name"
sha256sum "$archive" > "$archive.sha256"
echo "result_dir=$project_dir/$result_dir"
echo "result_archive=$project_dir/$archive"
echo "result_checksum=$project_dir/$archive.sha256"
