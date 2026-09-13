#!/usr/bin/env bash
# Short controller: audit completed raw files, submit only the next bounded array,
# then exit. A pending dependency consumes no CPUs.
#SBATCH --job-name=mj-n7-control
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

pipeline_cli=(julia --startup-file=no --project=. scripts/fast_ed_pipeline.jl)
current="$("${pipeline_cli[@]}" action --config="$BASE_CONFIG")"
IFS=$'\t' read -r current_action _ _ _ <<<"$current"
if [[ "$current_action" == "awaiting_results" ]]; then
    if ! decision="$("${pipeline_cli[@]}" advance --config="$BASE_CONFIG")"; then
        echo "Pipeline audit failed; state was changed to review." >&2
        exit 2
    fi
else
    decision="$current"
fi
IFS=$'\t' read -r action field1 field2 state_path <<<"$decision"
echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) pipeline action=$action state=$state_path"

case "$action" in
    solve)
        next_config="$field1"
        tasks="$field2"
        IFS=$'\t' read -r prepare_cpus prepare_threads configured_solve_cpus solve_threads configured_concurrency control_cpus \
            <<<"$("${pipeline_cli[@]}" resources --config="$BASE_CONFIG")"
        solve_cpus="${SOLVE_CPUS_OVERRIDE:-$configured_solve_cpus}"
        max_concurrent="${MAX_CONCURRENT_OVERRIDE:-$configured_concurrency}"
        [[ "$tasks" =~ ^[1-9][0-9]*$ && "$solve_cpus" =~ ^[1-9][0-9]*$ &&
           "$solve_threads" =~ ^[1-9][0-9]*$ && "$max_concurrent" =~ ^[1-9][0-9]*$ ]] || {
            "${pipeline_cli[@]}" review --config="$BASE_CONFIG" --reason=invalid_resource_fields
            exit 2
        }
        (( solve_cpus >= solve_threads )) || {
            "${pipeline_cli[@]}" review --config="$BASE_CONFIG" --reason=cpus_below_threads
            exit 2
        }
        array_last=$((tasks-1))
        if ! solve_raw=$(sbatch --parsable \
            --cpus-per-task="$solve_cpus" \
            --array="0-${array_last}%${max_concurrent}" \
            --export="ALL,PROJECT_ROOT=$PROJECT_ROOT,CONFIG=$next_config,THREADS=$solve_threads,FORCE_SOLVE=false" \
            slurm/fast_ed_sector_array.sbatch); then
            "${pipeline_cli[@]}" review --config="$BASE_CONFIG" --reason=solve_submission_failed
            exit 2
        fi
        solve_job="${solve_raw%%;*}"
        if ! controller_raw=$(sbatch --parsable \
            --dependency="afterany:$solve_job" \
            --cpus-per-task="$control_cpus" \
            --export="ALL,PROJECT_ROOT=$PROJECT_ROOT,BASE_CONFIG=$BASE_CONFIG,SOLVE_CPUS_OVERRIDE=$solve_cpus,MAX_CONCURRENT_OVERRIDE=$max_concurrent" \
            slurm/fast_ed_pipeline_controller.sbatch); then
            scancel "$solve_job"
            "${pipeline_cli[@]}" review --config="$BASE_CONFIG" --reason=controller_submission_failed
            exit 2
        fi
        controller_job="${controller_raw%%;*}"
        if ! "${pipeline_cli[@]}" record-submission --config="$BASE_CONFIG" \
            --stage="$(basename "$next_config" .toml)" \
            --solve-job="$solve_job" --controller-job="$controller_job" >/dev/null; then
            scancel "$controller_job" "$solve_job"
            "${pipeline_cli[@]}" review --config="$BASE_CONFIG" --reason=submission_record_failed
            exit 2
        fi
        echo "submitted solve=$solve_job controller=$controller_job tasks=$tasks max_concurrent=$max_concurrent cpus_per_task=$solve_cpus threads=$solve_threads"
        ;;
    bundle)
        archive="$field1"
        bundle_list="$field2"
        [[ "$archive" == "$PROJECT_ROOT"/* && "$bundle_list" == "$PROJECT_ROOT"/* ]] || {
            "${pipeline_cli[@]}" review --config="$BASE_CONFIG" --reason=bundle_path_outside_project
            exit 2
        }
        mkdir -p "$(dirname "$archive")"
        temporary="${archive}.tmp-${SLURM_JOB_ID:-$$}"
        trap 'rm -f -- "$temporary"' EXIT
        tar -C "$PROJECT_ROOT" -czf "$temporary" --files-from="$bundle_list"
        gzip -t "$temporary"
        mv -f -- "$temporary" "$archive"
        if command -v sha256sum >/dev/null 2>&1; then
            sha=$(sha256sum "$archive" | awk '{print $1}')
        else
            sha=$(shasum -a 256 "$archive" | awk '{print $1}')
        fi
        printf '%s  %s\n' "$sha" "$(basename "$archive")" >"${archive}.sha256"
        "${pipeline_cli[@]}" mark-bundled --config="$BASE_CONFIG" \
            --archive="$archive" --sha256="$sha" >/dev/null
        trap - EXIT
        echo "complete archive=$archive sha256=$sha"
        ;;
    review)
        echo "Pipeline stopped for review; no further ED work was submitted." >&2
        exit 3
        ;;
    complete)
        echo "Pipeline is already complete."
        ;;
    *)
        "${pipeline_cli[@]}" review --config="$BASE_CONFIG" --reason=unexpected_controller_action
        echo "Unexpected pipeline action: $action" >&2
        exit 2
        ;;
esac
