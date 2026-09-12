# N=7 fresh five-point scout, first new point

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

## Fresh restart: use the completed retained-N7 method

The active profile is `config/fast_ed/n7_uf0_165_k20_scout_restart.toml`:
`(Uf,Uf0,U0,Vf,Vf0,V0,t)=(0.46,1.65,4.14,0,0.55,0.34,0.5)`.
This is the first new N7 point; the completed retained center and N5/6 are unchanged.

The user explicitly requested a fresh start after job 548331 ran for over 3.5 hours.
Stop that job before rebuilding this point's cache. Do not reuse its matrix files or
sector results. `FORCE_REBUILD=true` rebuilds all four selected sector caches, and
`FORCE_SOLVE=true` forces fresh spectra. The new result run name is
`n7_uf0_165_scout_restart`. The retained center has a different cache identity and is
untouched. After this fresh preparation, all new scout tasks share the newly built
matrices: `H(mu)=H0+mu*Nf`. Rebuilding matrices separately for every mu is unnecessary.

Use the same staged method that completed the retained N7 point:

1. Prepare the current Hamiltonian's four sector caches once, using 8 CPUs/threads.
2. Run five scout mu values as 20 independent `(mu, sector)` Slurm array tasks:
   each task requests 8 CPUs/threads, BLAS 1; `--array=0-19%4` allows at most four
   simultaneous tasks. Each task exits and releases its own allocation when done.
3. Confirm all 20 tasks succeeded, then collect their small CSVs using one CPU.
   Do not attach collection to the original array's `afterok`: a failed index that
   is later rerun separately would leave that original dependency blocked.
4. Review the five-point q curve, actual spectra and competing branches. Only then
   choose the next small refinement grid, as with the old seven-point refinement.
   No automatic wide optimizer, automatic next Hamiltonian, or cache release.

The five scout values are
`[0.13512072184094, 0.14012072184094, 0.14512072184094, 0.15012072184094, 0.15512072184094]`.
Their center is `2*mu(N6)-mu(N5)`, using this point's audited N5/N6 values
`0.14128848931940/0.14320460558017`. N3/N4 guides are
`0.10853896038463/0.13537251915897`. These guides are not measured N7 muc values.

## Correctness and timing

The Hamiltonian, FuzzifiED, k=20, cold solver starts, five q terms and raw-rank
conventions are unchanged. Incomplete or identity-mismatched sector CSVs cannot be
collected as a valid full spectrum. The five-point minimum is only a scout candidate:
`collection_manifest.toml` explicitly records `stage="scout"`, `requires_review=true`,
and `muc_confirmed=false`. Even an interior sampled minimum is not a certified muc.
Inspect boundary behavior, all five residuals, factor, DeltaS/O, and state ordering;
choose refinement and any targeted rival-branch checks from those results. Do not
assign a high q to an unscored point or discard an inconvenient competing minimum.

The earlier retained-N7 preparation took 14m40s and its largest-sector pilot 7m59s.
Those measurements refer to one preparation and one sector, not an entire search.
The completed old workflow used five scout points followed by seven refinement
points. New-point queue time and total runtime have not been measured.

Commits `589656e` and `097cc0f` deviated from that workflow: the first held one
32-CPU allocation; the second serialized every sector inside one 8-CPU job. Both
also expanded the work to 9 local, 21 wide and 15 guard points plus multiple
refinements, with a 120-mu budget. The user rejected both deviations. The old
`scripts/fast_mu_search.jl` and `slurm/fast_mu_search.sbatch` entrypoints are removed.
`experimental/FastMuSearch.jl` and its profile/test remain historical prototypes,
not the active server workflow. Peak concurrency of four independent 8-CPU tasks
must not be confused with holding 32 CPUs for the lifetime of one long job.

## Server submission

First stop job 548331 and confirm it has left the queue. Pull `fast-ed-experiment`
in `/public/home/ruiqixu/MottJainED/MottJainED-fast-ed`, then run:

```bash
mkdir -p slurm-logs
export CONFIG=config/fast_ed/n7_uf0_165_k20_scout_restart.toml
export THREADS=8
PREP_RAW=$(sbatch --parsable --cpus-per-task=8 --export=ALL,FORCE_REBUILD=true slurm/fast_ed_prepare.sbatch)
PREP=${PREP_RAW%%;*}
printf '%s\n' "$PREP" > slurm-logs/last-n7-uf0165-prepare-job-id.txt
SCOUT_RAW=$(sbatch --parsable --cpus-per-task=8 --array=0-19%4 --dependency=afterok:"$PREP" --export=ALL,FORCE_SOLVE=true slurm/fast_ed_sector_array.sbatch)
SCOUT=${SCOUT_RAW%%;*}
printf '%s\n' "$SCOUT" > slurm-logs/last-n7-uf0165-scout-job-id.txt
squeue -j "$PREP,$SCOUT" -o '%.18i %.24j %.8T %.6C %.12M %R'
```

Only scout depends on the single preparation job; its tasks start after preparation
succeeds. Logs identify the config, task, mu/sector indices, threads and fresh-run
flags immediately; Julia logs dependency loading, operator preparation, cache writes,
matrix loading and the start of each eigensolve. Do not launch a second rebuild
while any task for this point is still running.

After all 20 array tasks are confirmed successful, collect explicitly:

```bash
sbatch --export=ALL,CONFIG=config/fast_ed/n7_uf0_165_k20_scout_restart.toml slurm/fast_ed_collect.sbatch
```

Download `output/fast_ed/runs/n7_uf0_165_scout_restart/`, plus preparation/array
logs and accounting statistics. It contains cache/settings-qualified
`scan_summary.csv`, `best_summary.csv`, `best_relations.csv`, `best_spectrum.csv`,
`collection_manifest.toml` and all per-mu small spectra. Do not download the large
matrix cache. Return these scout results for review before scheduling refinement.
The current profile forbids cache release; the retained cache and
`n7_k20_plot_files.tar.gz` remain preserved.

## Local validation

358 main-suite assertions and eight N3 integration assertions passed locally.
Slurm shell syntax and the mocked prepare/20-task launch checks also passed.

The main suite verifies cached/direct spectrum equivalence, forced cache and result
replacement, truncated-result rejection, the exact 20-task profile, and isolation
from the retained center. A separate N3 integration runs fresh preparation, all 20
sector solves and collection; it checks that the sampled minimum matches the N3
guide and that the collection remains explicitly preliminary:

```bash
julia --startup-file=no --project=. test/runtests.jl
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  julia --startup-file=no --project=. test/fast_ed_scout_restart_integration.jl
```

Shell launch checks use a mocked Julia executable to verify force flags, 8-thread
environment and the five-by-four task mapping. These are not real Slurm runtime or
performance measurements. No N7 ED is performed locally.
