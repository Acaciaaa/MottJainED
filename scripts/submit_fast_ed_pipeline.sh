#!/usr/bin/env bash

# Launch one bounded N=7 Hamiltonian pipeline. The script submits prepare, the
# five-point sector array, and a short dependent controller; it never waits while
# holding an allocation.
set -euo pipefail

if (( $# != 1 )); then
    echo "usage: bash scripts/submit_fast_ed_pipeline.sh CONFIG.toml" >&2
    exit 2
fi

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BASE_CONFIG="$1"
[[ "$BASE_CONFIG" = /* ]] || BASE_CONFIG="$PROJECT_ROOT/$BASE_CONFIG"
cd "$PROJECT_ROOT"
mkdir -p slurm-logs

pipeline_cli=(julia --startup-file=no --project=. scripts/fast_ed_pipeline.jl)
"${pipeline_cli[@]}" init --config="$BASE_CONFIG" >/dev/null
IFS=$'\t' read -r current_action initial_config initial_tasks state_path \
    <<<"$("${pipeline_cli[@]}" action --config="$BASE_CONFIG")"
[[ "$current_action" == "ready" ]] || {
    echo "pipeline is not ready (action=$current_action, state=$state_path)" >&2
    exit 2
}

IFS=$'\t' read -r configured_prepare_cpus prepare_threads configured_solve_cpus solve_threads configured_concurrency control_cpus \
    <<<"$("${pipeline_cli[@]}" resources --config="$BASE_CONFIG")"
prepare_cpus="${PREPARE_CPUS_OVERRIDE:-$configured_prepare_cpus}"
solve_cpus="${SOLVE_CPUS_OVERRIDE:-$configured_solve_cpus}"
max_concurrent="${MAX_CONCURRENT_OVERRIDE:-$configured_concurrency}"
for value in "$prepare_cpus" "$prepare_threads" "$solve_cpus" "$solve_threads" "$max_concurrent" "$control_cpus"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "invalid resource value: $value" >&2; exit 2; }
done
(( prepare_cpus >= prepare_threads && solve_cpus >= solve_threads )) || {
    echo "allocated CPUs must be at least the configured Julia threads" >&2
    exit 2
}

"${pipeline_cli[@]}" claim-launch --config="$BASE_CONFIG" >/dev/null
array_last=$((initial_tasks-1))
prepare_job=""
solve_job=""
controller_job=""

rollback_launch() {
    local job
    for job in "$controller_job" "$solve_job" "$prepare_job"; do
        [[ -z "$job" ]] || scancel "$job" >/dev/null 2>&1 || true
    done
    "${pipeline_cli[@]}" reset-launch --config="$BASE_CONFIG" >/dev/null 2>&1 || true
}

if ! prepare_raw=$(sbatch --parsable \
    --cpus-per-task="$prepare_cpus" \
    --export="ALL,PROJECT_ROOT=$PROJECT_ROOT,CONFIG=$initial_config,THREADS=$prepare_threads,FORCE_REBUILD=false" \
    slurm/fast_ed_prepare.sbatch); then
    rollback_launch
    exit 2
fi
prepare_job="${prepare_raw%%;*}"

if ! solve_raw=$(sbatch --parsable \
    --dependency="afterok:$prepare_job" \
    --kill-on-invalid-dep=yes \
    --cpus-per-task="$solve_cpus" \
    --array="0-${array_last}%${max_concurrent}" \
    --export="ALL,PROJECT_ROOT=$PROJECT_ROOT,CONFIG=$initial_config,THREADS=$solve_threads,FORCE_SOLVE=false" \
    slurm/fast_ed_sector_array.sbatch); then
    rollback_launch
    exit 2
fi
solve_job="${solve_raw%%;*}"

if ! controller_raw=$(sbatch --parsable \
    --dependency="afterany:$solve_job" \
    --cpus-per-task="$control_cpus" \
    --export="ALL,PROJECT_ROOT=$PROJECT_ROOT,BASE_CONFIG=$BASE_CONFIG,SOLVE_CPUS_OVERRIDE=$solve_cpus,MAX_CONCURRENT_OVERRIDE=$max_concurrent" \
    slurm/fast_ed_pipeline_controller.sbatch); then
    rollback_launch
    exit 2
fi
controller_job="${controller_raw%%;*}"

if ! "${pipeline_cli[@]}" record-submission --config="$BASE_CONFIG" \
    --stage=scout --prepare-job="$prepare_job" --solve-job="$solve_job" \
    --controller-job="$controller_job" >/dev/null; then
    rollback_launch
    exit 2
fi

printf 'prepare=%s\nscout=%s\ncontroller=%s\nstate=%s\n' \
    "$prepare_job" "$solve_job" "$controller_job" "$state_path"
printf 'resources: prepare=%sCPU/%sthread solve=%sCPU/%sthread max_concurrent=%s controller=%sCPU\n' \
    "$prepare_cpus" "$prepare_threads" "$solve_cpus" "$solve_threads" "$max_concurrent" "$control_cpus"
