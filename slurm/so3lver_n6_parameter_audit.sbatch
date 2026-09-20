#!/usr/bin/env bash
#SBATCH --job-name=so3-n6-audit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-logs/so3-n6-audit-%j.out

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
allocated_cpus="${SLURM_CPUS_PER_TASK:-1}"
export JULIA_NUM_THREADS="$allocated_cpus"
export OMP_NUM_THREADS="$allocated_cpus"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"

julia --startup-file=no --project=. scripts/so3lver_n6_parameter_audit.jl \
  --config=config/so3lver/n6_parameter_reference_audit.toml

result_parent="output/so3lver/parameter_search"
result_name="n6_reference_audit_v2"
archive="$result_parent/${result_name}_job-${SLURM_JOB_ID}.tar.gz"
tar -czf "$archive" -C "$result_parent" "$result_name"
sha256sum "$archive" > "$archive.sha256"
echo "result_archive=$project_dir/$archive"
echo "result_checksum=$project_dir/$archive.sha256"
