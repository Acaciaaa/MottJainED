#!/usr/bin/env bash
#SBATCH --job-name=so3-n6-conf-fix
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-logs/so3-n6-conf-fix-%j.out

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
result_name="n6_multistart_02_cg_corrected"
result_dir="$result_parent/$result_name"
mkdir -p "$result_dir"
cp "$config_path" "$result_dir/optimization_config.toml"

# Six deliberately separated four-dimensional starts.  Every worker builds an
# independent projected SO(3) workspace, evaluates four independently seeded
# Latin-hypercube points (24 total), runs one local optimizer, and writes only
# to its own signed checkpoint directory.
starts=(
  "0.48,1.68,0.3549772487682664,0.22881105175318384"
  "0.42,1.38,0.2199772487682664,0.16131105175318384"
  "0.54,1.98,0.4899772487682664,0.29631105175318384"
  "0.42,1.98,0.4899772487682664,0.29631105175318384"
  "0.54,1.38,0.4899772487682664,0.16131105175318384"
  "0.54,1.98,0.2199772487682664,0.16131105175318384"
)

pids=()
for lane in $(seq 1 6); do
  worker_dir="$result_dir/worker-$lane"
  worker_log="$result_dir/worker-$lane.stdout.log"
  mkdir -p "$worker_dir"
  echo "launching_worker=$lane start=${starts[$((lane - 1))]} log=$worker_log"
  srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 \
    julia --startup-file=no --project=. \
      scripts/so3lver_conformal_optimize.jl \
      --config="$config_path" \
      --output="$worker_dir" \
      --worker-index="$lane" \
      --parallel-workers=6 \
      --initial-samples=4 \
      --local-starts=1 \
      --forced-start="${starts[$((lane - 1))]}" \
    > "$worker_log" 2>&1 &
  pids+=("$!")
done

worker_failure=0
for index in "${!pids[@]}"; do
  lane=$((index + 1))
  if wait "${pids[$index]}"; then
    echo "worker_complete=$lane"
  else
    echo "worker_failed=$lane log=$result_dir/worker-$lane.stdout.log" >&2
    worker_failure=1
  fi
done
if (( worker_failure != 0 )); then
  echo "At least one conformal worker failed; retained all checkpoints in $result_dir" >&2
  exit 1
fi

julia --startup-file=no --project=. \
  scripts/so3lver_conformal_parallel_summary.jl \
  --config="$config_path" \
  --root="$result_dir" \
  --workers=6

archive="$result_parent/${result_name}_job-${SLURM_JOB_ID}.tar.gz"
tar -czf "$archive" -C "$result_parent" "$result_name"
sha256sum "$archive" > "$archive.sha256"
echo "result_dir=$project_dir/$result_dir"
echo "result_archive=$project_dir/$archive"
echo "result_checksum=$project_dir/$archive.sha256"
