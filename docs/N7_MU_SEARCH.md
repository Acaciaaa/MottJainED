# N=7 Uf0=1.65: two targeted follow-up points

## N=5/6 results audited on 2026-09-12

The downloaded `output/two_size_tuning/n56_retained_local_{uf0,vf0}_01/` directories
contain all 20 expected search rows and 786 individual q evaluations. Both CSV audits
pass; all 10 N=5/6 searches converged away from the boundary and their local/wide μ
differences are at most `1.73583283e-5`. Recomputing q from the five raw gaps agrees
to `4.44e-16`. The common center reproduces stage8 to less than `4e-11` in DeltaS/O.

| Varied parameter | Value | μc N5 | μc N6 | DeltaS N6−N5 | DeltaO N6−N5 | mean q |
|---|---:|---:|---:|---:|---:|---:|
| Uf0 | 1.65 | 0.14128849 | 0.14320461 | −0.02460426 | −0.13616659 | 0.12138697 |
| Uf0 | 1.834 | 0.14078291 | 0.14352012 | −0.00188866 | −0.18808935 | 0.13452760 |
| Uf0 | 2.00 | 0.14035232 | 0.14380869 | +0.02584737 | −0.22400969 | 0.18659178 |
| Vf0 | 0.45 | 0.14661204 | 0.15079875 | +0.01634306 | −0.21573760 | 0.16568066 |
| Vf0 | 0.65 | 0.13608279 | 0.13807345 | −0.01416549 | −0.15109460 | 0.13769402 |

These are local finite-size responses, not another Hamiltonian optimization. S drift
changes sign, while O drift stays negative. At Vf0=0.65, N6 also has a higher local
minimum at `mu=0.12771775, q=0.23693874`, separate from the best minimum at
`mu=0.13807345, q=0.11603076`. A single successful Brent run is insufficient evidence.

All five N6 searches had an unscored wide-Brent evaluation at `mu=0.0577864045`:
curlJ raw rank 3 was absent even at k=30. An unscored point cannot be treated as
evidence that q is large. Raw-rank conventions and the five score terms remain unchanged.

## Current result: refinement 548498 is complete, muc is still under review

All 28 tasks completed with exit code 0:0. The downloaded 28 CSV/TOML pairs contain
560 states with complete per-sector ranks 1–20. Cache/settings identities,
Hamiltonian, source hashes, solver settings, dimensions and eight-thread metadata
match the successful scout. Local Julia scoring and an independent Python raw-gap
reconstruction agree exactly at the saved precision.

| Fine-grid mu | q | DeltaS | DeltaO |
|---:|---:|---:|---:|
| 0.14425 | 0.09688635 | 1.69646420 | 2.82596127 |
| 0.14450 | 0.09698241 | 1.71441985 | 2.82346000 |
| 0.14475 | 0.09785799 | 1.73219906 | 2.82113748 |
| 0.14500 | 0.09945872 | 1.74978007 | 2.81898320 |
| 0.14525 | 0.10171605 | 1.76714455 | 2.81698725 |
| 0.14550 | 0.10455334 | 1.78427728 | 2.81514030 |
| 0.14575 | 0.10789152 | 1.80116575 | 2.81343359 |

The best sampled mu is the fine grid's left endpoint. The older scout has a higher
q at 0.14012072, but this leaves a large gap to the fine grid. Do not accept the
endpoint or move to another Hamiltonian yet. The five-point scout's quadratic
estimate near 0.1449 was too far right for the chosen fine interval.

Using the first three fine points, quadratic fits to q and q squared estimate
0.14434419 and 0.14434436. These estimates choose two new evaluations; neither is
a measured muc. The current provisional best has factor 0.02574342076446387;
its five scaled gaps are [1.04183428, 2.09844000, 3.07383779, 2.82770493, 3.01888586].

All seven ground states are singlets in (+,+). Fixed score ranks remain unchanged:
S raw rank 2, dS rank 1, J rank 1, curlJ raw rank 3, dJ rank 1 and T rank 1. J and
curlJ each have their expected pair of symmetry copies; switching the lowest copy
between Z/R sectors at machine precision is not evidence of a physical branch change.
The maximum required sector rank is 9. The largest quantum-number rounding error is
2.92e-10 and lowest-multiplet copy splitting is 4.89e-14. The factor and selected
gaps vary smoothly over the fine grid. Only one distinct O multiplet is visible;
these files cannot establish overlap continuity or exclude crossings with a higher O.

## Next calculation: only mu=0.14400 and 0.144375

Use `config/fast_ed/n7_uf0_165_k20_followup.toml`:

- 0.14400 checks whether q rises on the left of the current fine interval.
- 0.144375 samples close to the new estimated minimum, between existing 0.14425
  and 0.14450. Those existing neighbors give a 0.000125 local spacing.
- Two new mu values times four sectors means eight tasks, `--array=0-7%4`.
- Reuse the cache built by 548421, `output/fast_ed/cache/nm7_aaa8ccb4a8dd9fef`.
  No prepare, forced recalculation, dependency change or cache deletion is needed.
- Each independent task uses 8 CPUs/threads, BLAS 1; at most four run together.
  Hamiltonian, cold k=20 solver and five score terms are unchanged.
- Save results separately in `output/fast_ed/runs/n7_uf0_165_followup/`.

Review the new raw spectra **together with** the completed scout and seven-point
refinement. A two-point file alone cannot have an interior minimum. Check both
sides of the combined local minimum, fixed ranks, residuals, factor and S/O before
accepting an actual sampled mu. All profiles retain `allow_cache_release=false`.
No server collect job is requested; download the eight CSV/TOML pairs and logs.

The same-parameter N3/N4 guides and N5/N6 local/wide searches support the continuation
branch, but are not a global-minimum proof for N7. A boundary or conflicting branch
in the combined data calls for targeted review, not an automatic broad optimizer.
After this point passes review, the next planned Hamiltonian is Uf0=2.00 at
Vf0=0.55, V0=0.34. Later come Vf0=0.45 and 0.65 at Uf0=1.834. Generate and run one
Hamiltonian at a time; preserve the completed retained-center cache and results.

## Resources measured for 548498

| Measurement | Result |
|---|---:|
| Task elapsed range / median | 5m55s–8m44s / 6m31s |
| Largest MaxRSS | 15.28 GiB |
| First task start to estimated last task end | 53m35s |
| Allocated / actual CPU hours | 25.8800 / 8.6888 |
| CPU utilization of allocations | 33.57% |

The span uses log start timestamps plus sacct elapsed times; it excludes queue time
before the first task. For eight similar tasks with four concurrent slots, allow
roughly 12–20 minutes of execution plus queue/I/O delays; this is an estimate, not a
guarantee. Preserve the validated 8-CPU setup; no new thread-count benchmark was run.

Server Julia is 1.12.1 with FuzzifiED_jll 1.0.3+0; local scoring uses Julia 1.11.6
and a different JLL. Project/FuzzifiED matrix-source hashes match. Local analysis
reads the spectra only and preserves the server's recorded cache identity.

## Server commands: pull separately, then submit once

Before pull:

```bash
cd /public/home/ruiqixu/MottJainED/MottJainED-fast-ed
git switch fast-ed-experiment
mkdir -p slurm-logs
```

Pull separately; retry only this command on connection failures:

```bash
cd /public/home/ruiqixu/MottJainED/MottJainED-fast-ed
git pull --ff-only origin fast-ed-experiment
```

After successful pull, verify HEAD includes the commit reported with this change,
then submit only the follow-up array:

```bash
(
set -euo pipefail
cd /public/home/ruiqixu/MottJainED/MottJainED-fast-ed
test -s config/fast_ed/n7_uf0_165_k20_followup.toml
[ -z "$(squeue -h -u "$USER")" ] || { echo '队列中已有任务，未重复提交。'; exit 1; }
CACHE_DIR=output/fast_ed/cache/nm7_aaa8ccb4a8dd9fef
for part in zpos_rpos zpos_rneg zneg_rpos zneg_rneg; do
    test -s "$CACHE_DIR/$part.jld2" || { echo "缺少 cache：$part，停止提交。"; exit 1; }
done
mkdir -p slurm-logs
export PROJECT_ROOT="$PWD"
export CONFIG=config/fast_ed/n7_uf0_165_k20_followup.toml
export THREADS=8
N7_FOLLOWUP_RAW=$(sbatch --parsable --cpus-per-task=8 --array=0-7%4 --export=ALL,FORCE_SOLVE=false slurm/fast_ed_sector_array.sbatch)
N7_FOLLOWUP=${N7_FOLLOWUP_RAW%%;*}
printf '%s\n' "$N7_FOLLOWUP" | tee slurm-logs/last-n7-uf0165-followup-job-id.txt
squeue -j "$N7_FOLLOWUP" -o '%.18i %.24j %.8T %.6C %.12M %R'
)
```

After completion, return `output/fast_ed/runs/n7_uf0_165_followup/`, the eight
`mj-n7-sector-JOBID_*.out` logs, and sacct statistics. Keep the matrices on the server.

## Local audit artifacts

`output/fast_ed/analysis/n7_uf0_165_refine_548498/` contains the reproducible Julia
scoring script, independent Python check, resource/census check, seven-point and
combined twelve-point score tables, per-relation and selected-state tables, merged
spectra, raw-file SHA-256 hashes, `audit.toml`, resource statistics and original
`sacct.txt`. Its `best_sampled.csv` is provisional: `muc_confirmed=false` and
`requires_targeted_followup=true`. No N7 diagonalization was run locally.

The follow-up profile was loaded through FastED and checked against the completed
refinement: cache identity, solver settings, Hamiltonian and score agree; its two
mu values are new and its result directory is separate. All eight task mappings
were exercised through the existing Slurm script using a mock Julia executable;
they invoke only solve with eight threads and BLAS 1. Submission shell syntax and
the Git diff whitespace check pass. No solver or Slurm code changed.

## Retired workflow history

Job 548331 used a serial-sector automatic optimizer and was cancelled. Its cache and
small result directory were explicitly deleted before 548421 rebuilt the cache and
548422 reran the five-point scout. That fresh-start request does not mean rebuilding
the successful scout cache at each later mu refinement.

Commits 589656e/097cc0f deviated from the old small staged scans: fixed 32-CPU workers,
then fully serial 8-CPU execution, together with 9 local, 21 wide and 15 guard points
and a 120-mu budget. Those operational entrypoints were removed in bf834b6. The
historical FastMuSearch module/profile/test are not the active server path. The
current path is the existing independent sector array used for the retained N7
five-point scout and seven-point refinement.
