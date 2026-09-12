# N=7 cached chemical-potential search, first new point

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

## First point and search logic

Only `config/fast_ed/n7_uf0_165_search.toml` is prepared for a new N7 Hamiltonian:
`(Uf,Uf0,U0,Vf,Vf0,V0,t)=(0.46,1.65,4.14,0,0.55,0.34,0.5)`.
This is not the old stage11 compromise, whose Vf0/V0 differ.

Its same-parameter guide minima are N3 `0.10853896038463`, N4 `0.13537251915897`,
N5 `0.14128848931940`, N6 `0.14320460558017`. Use N6 as a search seed. The linear
extrapolation `0.14512072184094` is a guide, not a fixed chemical potential.

`experimental/FastMuSearch.jl` orchestrates the existing FSS scalar-search routines:

1. Local continuation around the N6 seed: half-width 0.02, 9 grid points, up to 3
   boundary expansions, valley refinement, and an independent whole-range Brent
   challenger over `[-0.095,0.305]`, as in the established workflow.
2. Independently evaluate a 21-point grid over that full range (spacing 0.02), then
   refine every discovered valid valley. This does not depend on the first result.
3. Independently evaluate a 15-point guard grid over `[0.095,0.165]` (spacing 0.005)
   and refine its valleys. This covers the guide minima, the predicted N7 neighborhood,
   and the old roughly 0.10 versus 0.14 branch ambiguity.
4. Select the lowest actually evaluated q across all searches. Check the winner's
   two sides at offsets 0.00025 and refine again. Mu tolerance is `2e-5`.
5. Repeat the selected spectrum with independent cold eigensolver starts, at the
   same k=20; require q and DeltaS/O differences below `1e-5`.

Every μ score requires complete, identity-matching sector CSVs. A truncated sector
file is recomputed. Shared μ evaluations use memoization and disk results across
phases/restarts. The 120-μ budget caps unexpected search growth; it is not a promise
that all 120 values will be needed. A final cold repeat is one additional solve.

Disagreement, a missing score inside the dense guard, a winner outside the guard,
failed convergence, or a failed cold repeat leaves the data/cache and exits with a
review-needed error. A successful local/grid audit is **not** a mathematical proof
of a global minimum: unsampled narrow valleys and unscored outer regions remain
explicit limitations. The audit always records `global_minimum_proven=false`.
Do not turn a finite-resolution score minimum into a claim of a thermodynamic critical point.

## Cache and resources

Use `H(mu)=H0+mu*Nf` with the existing per-sector H0/Nf/L2/C2 JLD2 cache. Build this
one Hamiltonian's cache once (approximately 25 GB at N7). Four local Julia processes
each load one sector once and retain its matrices while the controller selects μ.
Each worker uses 8 threads, BLAS uses 1, and the Slurm job requests 32 CPUs on one
node. The controller waits while the four sectors are solved in parallel. This is
the same four-sector parallelism as the earlier array, now with matrices resident
across μ values. It has not yet been resource-profiled at N7 as one combined job;
inspect its MaxRSS/CPU usage after the first run.

Cold eigensolver starts are intentional: reusing an old eigenvector is not part of
the speed claim. The gains come from matrix reuse, persistent processes, sector
parallelism, and reusing completed μ results. FuzzifiED and the Hamiltonian are unchanged.

Preparation, search and collection run in one Slurm job. No later Hamiltonian is
submitted, and no cache is deleted automatically. The already-completed central N7
cache and `n7_k20_plot_files.tar.gz` are preserved.

## Server entry and results

From `/public/home/ruiqixu/MottJainED/MottJainED-fast-ed`, after pulling
`fast-ed-experiment` and creating `slurm-logs`:

```bash
sbatch slurm/fast_mu_search.sbatch
```

The log prints the exact cache/settings-qualified result directory under
`output/fast_ed/runs/n7_uf0_165_search/`. Download that small run directory, not
`output/fast_ed/cache/`. It includes:

- `mu_search_evaluations.csv`: full μ/q trace, with invalid scores and reasons;
- `mu_search_audit.toml`: coverage, convergence, competing-branch comparisons and issues;
- `scan_summary.csv`, `best_summary.csv`, `best_relations.csv`, `best_spectrum.csv`;
- `cold_repeat.csv`, `collection_manifest.toml`, `evaluated_profile.toml`;
- per-μ small sector and merged-spectrum CSVs for spectral continuity/rank checks.

The static input profile's `fast_ed.mus` contains only the initial seed; the search
chooses additional values dynamically. For later collection use the saved
`evaluated_profile.toml`, which lists the actual evaluated values. Do not collect
with the original one-seed profile and overwrite the full summaries.

`FastED.validate_collection` still rejects any invalid q in a full scan. A
`mu_search_audit.toml` that passes local/grid checks while listing unscored outer
points does not override that restriction or authorize cache release. Analyze the
downloaded spectrum/trace first; only after resolving coverage limits should an
explicit accepted-result/cache-release step be prepared. Failures retain the cache.

Re-submit the same script after a stopped/failed job to replay the search and reuse
completed sector outputs. Do not run two copies of this same point concurrently.

## Validation completed before the first server run

On local Julia 1.11.6, all 346 main-suite tests and 12 bounded integration checks pass:

```bash
julia --startup-file=no --project=. test/runtests.jl
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  julia --startup-file=no --project=. test/fast_mu_search_integration.jl
```

The integration test constructs only temporary N3 matrices, starts four one-thread
workers, runs the full search/collection/cold-repeat path, and replays it to verify
disk-result reuse. Its minimum `0.10854498882565` agrees with the downloaded N3 guide
within `6.03e-6`. Unit checks also exercise a hidden competing valley near 0.10,
an unscored point inside the guard, boundary and budget failures, resident/direct
spectrum equivalence, and truncated-CSV recovery. Slurm syntax and the 32-CPU guard
were checked separately. No new N7 ED calculation was run locally.
