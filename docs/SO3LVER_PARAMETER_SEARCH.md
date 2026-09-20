# SO(3)lver N=6 parameter search

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
output/so3lver/parameter_search/n6_reference_audit_v2/audit_summary.toml
output/so3lver/parameter_search/n6_reference_audit_v2/scores.csv
output/so3lver/parameter_search/n6_reference_audit_v2/state_tracking.csv
output/so3lver/parameter_search/n6_reference_audit_v2/normalized_jacobian.csv
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

The revised audit therefore keeps all six generator competitors and raises the
one-time conventional calculation to `k=80`. Its local `V0` difference is
reduced to 0.02 so that the Jacobian tests the anchor branch instead of
straddling the rearrangement. The Slurm wrapper packages the four result files
into `n6_reference_audit_v2_job-<jobid>.tar.gz` with a SHA-256 sidecar; download
that archive rather than pasting the files into a terminal transcript.

## Gated N=6 optimization

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
`output/so3lver/parameter_search/n6_multistart_02/`. The Slurm wrapper also
creates `n6_multistart_02_job-<jobid>.tar.gz` and its SHA-256 sidecar. In particular,
`best.toml` is only a finite-N candidate.  It must pass the unused `boxO/boxJ`,
generator, density/phase, N=7, and ultimately N=8 holdouts before it can replace
the baseline.
