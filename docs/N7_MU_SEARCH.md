# N=7 Uf0=1.65: completed scout and seven-point refinement

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

## Current step: seven-point refinement (2026-09-13)

The user downloaded the raw spectra and requested local analysis without a server
collection job. Preparation 548421 and all 20 tasks of scout array 548422 completed
with exit code 0:0. The cache manifest and 20 CSV/TOML pairs identify the intended
Hamiltonian `(Uf,Uf0,U0,Vf,Vf0,V0,t)=(0.46,1.65,4.14,0,0.55,0.34,0.5)`.
This is not the stage11 compromise with different Vf0/V0.

| Scout mu | q | DeltaS | DeltaO |
|---:|---:|---:|---:|
| 0.13512072184094 | 0.23075022 | 1.57561738 | 3.16332619 |
| 0.14012072184094 | 0.18052423 | 1.41758253 | 2.90544710 |
| 0.14512072184094 | 0.10047138 | 1.75819320 | 2.81800020 |
| 0.15012072184094 | 0.19312703 | 2.05289253 | 2.79942730 |
| 0.15512072184094 | 0.27789490 | 2.25381571 | 2.80376954 |

All five points contain 80 states and all five q relations. All 20 CSV/TOML pairs
pass identity, complete-rank, solver-setting and dimension checks; independent
Python reconstruction agrees with the Julia q/factor/DeltaS/O values exactly at the
saved precision. Maximum quantum-number rounding error is 2.50e-10, and the largest
lowest-multiplet copy splitting is 4.80e-14. Scoring uses the existing
raw-rank convention: S raw rank 2, J raw rank 1, curlJ raw rank 3 (the second distinct
multiplet), dJ raw rank 1, and T raw rank 1. No deduplication precedes this scoring.
The largest required sector rank over the scout is 13; at its best sampled point it
is 9. Cold k=20 starts, eigensolver tolerances and the five score terms stay unchanged.

The sampled minimum is interior. Fits through its two neighboring points estimate
mu=0.14493829 from q and 0.14488384 from q squared. These are only guides for choosing
the next grid. The existing N3/N4 full-grid and N5/N6 local/wide checks at this same
Hamiltonian select the continuation branch leading toward this neighborhood; they do
not prove that N7 has no unsampled competing branch. In particular, an old minimum
near 0.10 at different couplings must not be transplanted to this point.

Use `config/fast_ed/n7_uf0_165_k20_refine.toml` to compute exactly:

```text
0.14425, 0.14450, 0.14475, 0.14500, 0.14525, 0.14550, 0.14575
```

The spacing is 0.00025, as in the completed retained-N7 refinement. All seven values
are new. After receiving these results, compare the actual minimum and both sides,
five residuals, factor, scalar ordering, lowest O and the scout/guide branches. A
boundary minimum, discontinuity, missing level or conflicting branch calls for
specific additional points; do not automatically launch a broad optimizer. Do not
report an interpolated estimate as a measured muc or a proven critical point.
The raw scout includes only the lowest O multiplet's two copies, so it does not
establish continuity of a second, higher O level.

## Resources: measured scout and unchanged next-run settings

| Measurement | Result |
|---|---:|
| Preparation elapsed / MaxRSS | 13m51s / 22.99 GiB |
| Sector task elapsed range / median | 5m34s–9m28s / 6m50s |
| Largest sector-task MaxRSS | 14.13 GiB |
| Scout span, first task start to last task end | about 38m43s |
| Preparation start to last scout task end | about 52m34s |
| Scout allocated CPU hours / actual CPU hours | 18.6822 / 6.1280 |

Spans use the log start timestamps and sacct elapsed times; they exclude queue time
before preparation and are not a guaranteed runtime for refinement. Scout CPU
utilization averaged 32.80% of the 8-CPU task allocations. No new thread-count tuning
has been measured, so preserve the requested and previously validated 8-thread setup.

Keep the matrices freshly built by 548421. The server cache is
`output/fast_ed/cache/nm7_aaa8ccb4a8dd9fef`; refinement's unchanged Hamiltonian,
source and solver definitions retain the same cache/settings identities. Do not
prepare again, delete this cache, or use FORCE_SOLVE=true for normal continuation.
The 28 independent tasks each request 8 CPUs/threads, BLAS 1, with at most four
simultaneous tasks. Each exits and releases its own allocation. The retained center
cache `nm7_20d7825bc607cd05` stays preserved, as do the downloaded scout results.

The server manifest records Julia 1.12.1 and FuzzifiED_jll 1.0.3+0, while local
analysis uses Julia 1.11.6 with a different JLL. Matrix-defining project/FuzzifiED
source hashes match. Local work reads CSVs and computes q only; it does not generate
N7 eigenvalues or replace the server's recorded cache identity with the local one.
Do not update dependencies during this refinement.

## Server commands: keep pull separate from submission

Before pull, only select the saved branch and ensure the log directory exists:

```bash
cd /public/home/ruiqixu/MottJainED/MottJainED-fast-ed
git switch fast-ed-experiment
mkdir -p slurm-logs
```

Pull separately; repeat only this command on connection failures:

```bash
git pull --ff-only origin fast-ed-experiment
```

After a successful pull, submit only the refinement array (no prepare, dependency,
or server collection job). Check that HEAD includes the refinement commit reported
with this change before submitting:

```bash
export PROJECT_ROOT="$PWD"
export CONFIG=config/fast_ed/n7_uf0_165_k20_refine.toml
export THREADS=8
N7_REFINE_RAW=$(sbatch --parsable --cpus-per-task=8 --array=0-27%4 --export=ALL,FORCE_SOLVE=false slurm/fast_ed_sector_array.sbatch)
N7_REFINE=${N7_REFINE_RAW%%;*}
printf '%s\n' "$N7_REFINE" | tee slurm-logs/last-n7-uf0165-refine-job-id.txt
squeue -j "$N7_REFINE" -o '%.18i %.24j %.8T %.6C %.12M %R'
```

After completion, return the entire small directory
`output/fast_ed/runs/n7_uf0_165_refine/`, its `mj-n7-sector-JOBID_*.out` logs and
sacct resource statistics. Keep the cache on the server for any needed follow-up.
No server summary is required; analysis is done locally from the raw CSV/TOML files.

## Audit artifacts and validation

Local derived files are under `output/fast_ed/analysis/n7_uf0_165_scout_548422/`:
`scan_summary.csv`, `relations.csv`, `selected_states.csv`, per-mu merged spectra,
`audit.toml`, `source_hashes.csv`, `resource_audit.json`, `resources.csv`, and the
original pasted `sacct.txt`. Raw downloaded files are preserved. The recorded source
hashes, complete per-sector ranks, metadata and independent five-gap reconstruction
allow the numerical conclusions to be checked without transferring large matrices.

All 372 main-suite assertions pass. Profile regression tests check the unchanged
Hamiltonian/cache/settings/score,
seven new mu values at 0.00025 spacing, 28 tasks, an isolated result directory and
review-before-acceptance/cache-release flags. Shell checks exercise all 28 task
indices through the existing array script. No N7 ED is run locally.

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
