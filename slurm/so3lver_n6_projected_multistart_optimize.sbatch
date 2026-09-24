#!/usr/bin/env bash
#SBATCH --job-name=so3-n6-proj-par
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=sdicnormal
#SBATCH --output=slurm-logs/so3-n6-proj-par-%j.out

set -euo pipefail
project_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$project_dir"

result_parent="output/so3lver/parameter_search"
result_name="n6_projected_parallel_multistart_01"
result_dir="$result_parent/$result_name"
config_path="config/so3lver/n6_projected_multistart_search.toml"
mkdir -p "$result_dir"
cp "$config_path" "$result_dir/search_config.toml"

# Six deliberately separated first-round starts.  Each worker also evaluates
# five independently seeded Latin-hypercube points, for 30 exploratory points
# in total.  Workers never write to the same checkpoint or CSV.
starts=(
  "0.48,1.68,0.3549772487682664"
  "0.42,1.38,0.2199772487682664"
  "0.54,1.98,0.4899772487682664"
  "0.42,1.98,0.4899772487682664"
  "0.54,1.38,0.4899772487682664"
  "0.54,1.98,0.2199772487682664"
)

pids=()
for lane in $(seq 1 6); do
  worker_dir="$result_dir/worker-$lane"
  worker_log="$result_dir/worker-$lane.stdout.log"
  mkdir -p "$worker_dir"
  echo "launching_worker=$lane start=${starts[$((lane - 1))]} log=$worker_log"
  srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 \
    julia --startup-file=no --project=. \
      scripts/so3lver_n6_projected_multistart_optimize.jl \
      --config="$config_path" \
      --output="$worker_dir" \
      --worker-index="$lane" \
      --parallel-workers=6 \
      --initial-samples=5 \
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
  echo "At least one parallel worker failed; retained all checkpoints in $result_dir" >&2
  exit 1
fi

julia --startup-file=no --project=. \
  scripts/so3lver_projected_parallel_summary.jl \
  --config="$config_path" \
  --root="$result_dir" \
  --workers=6

archive="$result_parent/${result_name}_job-${SLURM_JOB_ID}.tar.gz"
tar -czf "$archive" -C "$result_parent" "$result_name"
sha256sum "$archive" > "$archive.sha256"
echo "result_archive=$project_dir/$archive"
echo "result_checksum=$project_dir/$archive.sha256"
