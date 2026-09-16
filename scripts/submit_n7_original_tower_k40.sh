#!/usr/bin/env bash

# Submit the bounded N=7 k=40 standard-family tower workflow.  This deliberately
# has no prepare job: it requires and reuses the completed original-Hamiltonian
# N=7 matrix cache already present on the server.
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONFIG="${CONFIG:-config/fast_ed/n7_original_tower_k40.toml}"
SERVER_CACHE="${SERVER_CACHE:-output/fast_ed/cache/nm7_6ff20bdae301b353}"
cd "$PROJECT_ROOT"
mkdir -p slurm-logs

[[ -f "$SERVER_CACHE/cache_manifest.toml" ]] || {
    echo "required completed cache is missing: $SERVER_CACHE" >&2
    exit 2
}
grep -Eq '^cache_id = "6ff20bdae301b353"$' "$SERVER_CACHE/cache_manifest.toml" || {
    echo "cache manifest has the wrong cache_id: $SERVER_CACHE" >&2
    exit 2
}
grep -Eq '^complete = true$' "$SERVER_CACHE/cache_manifest.toml" || {
    echo "cache manifest is not complete: $SERVER_CACHE" >&2
    exit 2
}

array_raw=$(sbatch --parsable \
    --array=0-3%4 \
    --export="ALL,PROJECT_ROOT=$PROJECT_ROOT,CONFIG=$CONFIG,THREADS=8,FORCE_SOLVE=false" \
    slurm/fast_ed_vector_sector_array.sbatch)
array_job="${array_raw%%;*}"

collect_raw=$(sbatch --parsable \
    --dependency="afterok:$array_job" \
    --kill-on-invalid-dep=yes \
    --export="ALL,PROJECT_ROOT=$PROJECT_ROOT,CONFIG=$CONFIG,THREADS=8" \
    slurm/fast_ed_tower_collect.sbatch)
collect_job="${collect_raw%%;*}"

printf 'vector_array=%s\ncollector=%s\n' "$array_job" "$collect_job"
printf 'cache=%s\nconfig=%s\n' "$SERVER_CACHE" "$CONFIG"
printf 'results=output/fast_ed/generator/n7_original_muc_k40\n'
