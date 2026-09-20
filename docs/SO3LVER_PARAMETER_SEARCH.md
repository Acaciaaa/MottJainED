# SO(3)lver N=6 parameter search

## Current recommended search: nested stable-six objective

The active N=6 search no longer uses `boxS` or a microscopic-generator gate.
The audit showed that `ddS` remains rank- and overlap-stable throughout the
local probes, whereas `boxS` can exchange rank across the useful `V0` range.
The current training score therefore contains exactly

```text
dS-S, ddS-dS, J, curlJ, dJ(rank1), T(rank1).
```

`Uf0`, `Vf0`, and `V0` are optimized jointly.  `mu` is not optimized once and
then held fixed: every outer point receives a complete coarse scan over the
configured `mu` interval, followed by independent Brent refinement of every
valid local valley (up to the configured limit).  This is the key guard
against mistaking a sequentially tuned point for a joint optimum.

The outer search first uses deterministic Latin-hypercube coverage, then runs
several diverse three-dimensional Nelder--Mead starts.  If the best point lies
against the current trust-region boundary, the region is recentered and
expanded for one more round.  Every scored state must keep its fixed physical
rank and pass the direct overlap gate.  Line continuation from the anchor is
allowed only as an identity fallback; it never relabels a rank.  A final full
`mu` reprofile and six one-sided outer-parameter neighbors are included in the
same job, so there is no separate generator or audit job to launch.

Run on the server with:

```bash
mkdir -p slurm-logs
sbatch slurm/so3lver_n6_nested_optimize.sbatch
```

The wrapper requests eight CPUs on `sdicnormal` and deliberately specifies
neither memory nor a time limit.  It packages the complete result directory as
`n6_nested_six_term_01_job-<jobid>.tar.gz`, with a SHA-256 sidecar.  Download
that archive for analysis instead of pasting long CSV or log output.

The main decision file is `best.toml`.  `accepted = true` requires an interior
`mu`, preservation of the original five-term locator to the configured ratio,
and no improvement at any of the six local robustness probes.  Regardless of
that flag, `top_candidates.csv`, `robustness_neighbors.csv`, and the complete
point/residual/identity traces are retained.  The result is still only an N=6
candidate; N=7 and N=8 are later validation stages, not part of this objective.

## Historical seven-term audit and optimizer

The multi-parameter search is deliberately split into a reference audit and a
later optimization.  Do not enable the provisional seven-relation score merely
because the required eigenvalues are present.

## Why the audit is mandatory

The fixed five-relation score fits one common scale from five gaps, leaving only
four independent dimensionless conditions.  Optimizing
`Uf0,Vf0,V0,mu` against those conditions at one size can therefore overfit the
finite-size spectrum.  Two scalar relations are useful only if their states do
not silently exchange rank:

```text
ddS-dS   target 1   (singlet L=2 physical rank 2 minus singlet L=1 rank 1)
boxS-S   target 2   (singlet L=0 physical rank 3 minus rank 2)
```

The reference audit provides three independent gates:

1. At the central N=6 point, a conventional Fock-basis calculation fits the
   microscopic generator only from `S -> dS`.  Applying the frozen generator
   back to the scalar tower checks whether `boxS` and `ddS` are the leading
   expected non-parent components and whether the two expected low-energy
   subspaces capture the generated vectors.
2. The SO(3)lver solver evaluates the center and positive/negative perturbations
   of all four proposed free parameters.  Direct eigenvector overlaps must map
   every scored state back to the same physical rank.  A rank exchange or an
   overlap below the configured threshold fails the audit.
3. The normalized finite-difference residual Jacobian must have rank four.  Its
   singular values and condition number are reported so correlated standards
   are visible rather than hidden by their count.  The checked-in profile also
   rejects a normalized condition number above 1000.

No optimization is launched by this audit.  The separate optimizer refuses to
start unless `audit_summary.toml` reports `passed = true` and its configuration
and source hashes still match.

## Command

From the `so3lver-experiment` worktree:

```bash
julia -t 8 --project=. scripts/so3lver_n6_parameter_audit.jl \
  --config=config/so3lver/n6_parameter_reference_audit.toml
```

The checked-in profile fixes the original audited baseline and its accepted
N=6 chemical potential.  It writes:

```text
output/so3lver/parameter_search/n6_reference_audit_v3/audit_summary.toml
output/so3lver/parameter_search/n6_reference_audit_v3/scores.csv
output/so3lver/parameter_search/n6_reference_audit_v3/state_tracking.csv
output/so3lver/parameter_search/n6_reference_audit_v3/normalized_jacobian.csv
```

The process still exits normally when a scientific gate fails, but prints
`AUDIT_NOT_CERTIFIED` and records `passed = false`.  This keeps a negative
scientific result distinct from a scheduler or numerical crash.  Inspect
`audit_summary.toml` instead of weakening the thresholds or changing ranks to
force a pass.

The first server audit on 2026-09-20 completed normally but was not certified.
It exactly reproduced the five-relation anchor score and found a full-rank
seven-residual Jacobian, but the conventional `k=20` spectrum contained only
four of the six requested singlet L=0 levels. More importantly, the original
`V0` probe at `0.525-0.08=0.445` crossed a low-scalar rearrangement: the anchor
ground-state overlap fell to 0.533, `S` mapped from rank 2 to rank 1, and
`boxS` mapped from rank 3 to rank 4. This is a real branch-safety failure, not
a reason to relabel the states.

The second audit kept all six requested generator competitors and raised the
one-time conventional calculation to `k=80`. Its local `V0` difference is
0.02, so the Jacobian tested the anchor branch instead of straddling the
rearrangement. All eight parameter probes then preserved every fixed rank; the
smallest expected overlap was 0.899. The residual Jacobian had rank 4/4,
singular values `(0.4086,0.03628,0.01301,0.002921)`, and condition number
139.9. The generator still stopped before overlap analysis because the `k=80`
spectrum contained five rather than six L=0 singlets.

The current gate requests up to six resolved competitors but no longer treats
an arbitrary exact count as physics. It requires at least four L=0 and three
L=2 resolved levels. All omitted states are represented by the unresolved
overlap `1-sum(resolved overlaps)`, which is used as a single worst-case
competitor. Thus `boxS` and `ddS` pass only when their individual overlaps beat
even the total omitted weight, in addition to the existing absolute-subspace,
resolved-fraction, and fit-fidelity thresholds. This is stricter than simply
dropping the sixth state. The Slurm wrapper packages the result as
`n6_reference_audit_v3_job-<jobid>.tar.gz` with a SHA-256 sidecar; download that
archive rather than pasting files into a terminal transcript.

## Historical gated N=6 optimization

Only after the audit passes, run:

```bash
julia -t 8 --project=. scripts/so3lver_n6_parameter_optimize.jl \
  --config=config/so3lver/n6_parameter_reference_audit.toml
```

The first profile uses three local, identity-gated Nelder--Mead starts for
`Uf0,Vf0,V0,mu`, while `Uf=.46`, `U0=4.14`, `Vf=0`, and `t=.5` remain fixed.
The simplex coordinates are scaled separately in physical units for each
parameter. In particular, the broad allowed `V0` interval no longer turns a
single generic normalized step into a `V0=0.09` jump across the scalar
rearrangement.
Every evaluated point is rejected if a scored eigenvector no longer maps to the
same anchor rank with the configured minimum overlap, or if the global ground
branch is no longer singlet L=0.  The scalar objective is the seven-relation RMS
plus a configured penalty on the largest individual residual.  A final wide
mu grid is a branch guard; disagreement is recorded instead of silently
accepting a local mu valley.

The optimizer writes its full checkpoint trace under
`output/so3lver/parameter_search/n6_multistart_03/`. The Slurm wrapper also
creates `n6_multistart_03_job-<jobid>.tar.gz` and its SHA-256 sidecar. In particular,
`best.toml` is only a finite-N candidate.  It must pass the unused `boxO/boxJ`,
generator, density/phase, N=7, and ultimately N=8 holdouts before it can replace
the baseline.
