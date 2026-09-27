# Projected conformal-algebra optimization

## Current implementation: full-algebra fit (2026-09-27)

The production profile now writes `n6_multistart_03_full_algebra`. It is not
compatible with either previous search trace. The CG norm correction remains
in place; the subsequent audit identified and corrected different problems:

1. The eight old pair-projector vector operators do not span all local
   Hamiltonian-density moments. Seven actual local-density moments have been
   added, including the light-Laplacian/heavy contact and heavy density.
2. The generalized eigenproblem fitted a proxy while the outer search scored
   the full commutators. It is now only an initializer: coefficients, their
   common amplitude, and the cylinder scale are jointly refined against
   exactly the objective used by the Hamiltonian optimizer.
3. Conserved-current shortening is tested through the `L -> L-1` component
   of `P`: on `J` in training and on `T` in holdout. Generic primary algebra
   alone does not require the protected dimensions of these currents.
4. A global overlap cutoff against a remote anchor is replaced by adaptive
   same-rank continuation. Smooth deformation is allowed; detected rank
   crossings and ground-branch changes in computed sectors are rejected.
   Global anchor overlap is retained only as a diagnostic.

### Generator, representation, and projection

The calculation continues to use only native `FuzzifiED.SO3lver` composite
spaces and direct Jack/Laughlin-1/3 projection, not the old full-Fock solver.
A single coefficient vector is shared by singlet and adjoint representations.
The `local_density` basis retains the old eight candidates and appends
`density_Uf, density_Vf, density_U0, density_V0, density_Uf0, density_Vf0,
density_heavy`. These are actual rank-one moments of the local density and
Laplacian products used in the scalar Hamiltonian, plus heavy density.
Projection-induced redundancies are removed using coupling-independent
full-block probes (15 candidates, numerical rank 10 in the N4 regression).
This is not an exhaustive basis of every possible local improvement.

The definitions remain `D=(H-E0)/factor`, `Lambda=P+K`,
`P=(Lambda+[D,Lambda])/2`, `K=(Lambda-[D,Lambda])/2`, with physical reduced-norm
weight `(2L'+1)/(2L+1)`. All intermediate states in two-operator products are
in complete projected angular-momentum blocks, not a low-energy truncation.
The exact full-energy definition can amplify UV transitions; sensitivity to
that choice is still a separate scientific question, not resolved by this fix.

### Identical inner and outer objectives

The training sources are `S,O,J,dS,curlJ`. `dS` is singlet `L=1`, rank 1;
`curlJ` is adjoint `L=1`, rank 2. These names are fixed-rank identification
hypotheses, not a proof of their continuum identity. They test algebra closure,
not `K|descendant>=0`. `T` is excluded from the initializer, full inner fit,
and outer objective. Its state identity is monitored and its algebra is an
independent final holdout veto.

Each training primary contributes `K-primary`, dilation, mixed `[K,P]`,
`[P,P]`, `[K,K]`, and low-energy leakage. Each descendant contributes these
except `K-primary`; `J` adds shortening, and the vacuum contributes one
constraint. The objective is the weighted mean of these 30 squared-norm
ratios plus `worst_weight` times their largest value. Leakage has weight 0.25;
the other terms have weight 1. There is no integer-gap or spectrum-score
term. The outer coordinates remain `Uf,Uf0,Vf0,muc`, with `U0=9Uf`,
`V0=Vf=0`, `t=0.5`.

The full inner fit uses small Gram matrices precomputed from complete-block
operator actions, then jointly varies direction, amplitude and scale with
two deterministic local starts. At every Hamiltonian evaluation the compact
objective is checked against independent direct operator application.
This is a local nonlinear fit, not a guaranteed global coefficient minimum.
Nonconverged fits are recorded and cannot receive final search acceptance.

Same-rank continuation uses straight-line anchor-to-point paths, at least
0.90 squared overlap per accepted step, and step sizes at most 0.10 of the
hard-bound widths. Failed steps are bisected (depth 10, up to 64 extra spectra
per point). Ground ordering is checked only in computed sectors. Crossings
between sampled steps or in uncomputed sectors cannot be ruled out.

### N6 preflight: one CPU, no parameter search

```bash
mkdir -p slurm-logs
sbatch slurm/so3lver_n6_conformal_pilot.sbatch
```

This evaluates the anchor and the previous corrected-search winner, then
cold-repeats the latter. It does not submit any follow-up job. Inspect
`pilot.toml`, the evaluation trace (per-point time and inner convergence),
and residual/coefficient tables before launching the full search.
`pilot_passed` means numerical convergence and repeatability, not a good CFT
or passing state identity at the old winner. Peak RSS is also recorded.
Download both files from `output/so3lver/conformal_optimization/`:

- `n6_full_algebra_pilot_job-JOBID.tar.gz`
- `n6_full_algebra_pilot_job-JOBID.tar.gz.sha256`

Archives are written even on ordinary pilot failure, so diagnostics survive.
N6 runtime, convergence, and physical results remain unverified locally.

### Full search: six independent one-CPU workers

```bash
sbatch slurm/so3lver_n6_conformal_optimize.sbatch
```

The job still runs six Julia processes simultaneously, with six deliberately
separated starts and 24 total Latin-hypercube points per exploration round.
If a prescribed start fails a gate, the worker first backs off along its own
direction instead of silently duplicating an anchor start. Each worker has
separate signed checkpoints. Neither Slurm file requests `--mem` or `--time`:
CPU count determines memory on this cluster. Measure the expanded-basis cost
with the pilot; do not assume the old implementation's roughly 2 GiB usage.
No projected-spectrum jobs or their files are changed by either entry point.

Final acceptance requires cold repeatability, valid identities and scale,
inner-fit convergence, coordinate-neighbor checks, independent-start consensus,
hard-bound clearance, and the `T` holdout veto. Legacy spectral scores are
post-selection diagnostics only. `algebra_fit_loss` is retained as the legacy
initializer diagnostic, not the objective being optimized; use
`inner_objective_initial/final`, `inner_converged`, and `objective` instead.
Best-point outputs now include actual generator coefficients and channel
residuals. Source hashes cover the driver, operator module, full-fit module,
tracking module, and configuration, and are captured at process startup.

Download the full-search archive and checksum:
`n6_multistart_03_full_algebra_job-JOBID.tar.gz` and `.tar.gz.sha256` from
the same parent directory. Neither prior result directory is overwritten.

### Local regression (N4 only)

```bash
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project=. --threads=1 test/so3lver_runtests.jl
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project=. --threads=1 test/so3lver_conformal_fit_runtests.jl
```

Tests cover exact free-scalar algebra including spin, projected generator
Hermiticity, the previously missing density moment, continuation versus rank
crossings, holdout exclusion, and compact/direct objective agreement.
Passing these and a small search test establishes implementation consistency,
not that an interacting CFT or better N6 point has been found.

---

## Historical v1/v2 notes — not current optimizer defaults

The remainder preserves earlier method descriptions and numerical results.
In particular, references below to an eight-candidate fit, a three-source
training set, or a global anchor-overlap gate describe old implementations.
The two historical fixed-point diagnostic scripts still use their legacy
profiles for reproducibility; use the N6 preflight above for the new method.

> **September 2026 normalization correction.**  An independent operator audit
> verified that the native SO(3)lver rank-one tensors, their Hermitian phases,
> and the direct Jack/Laughlin projection are correct.  It found a later
> reduced-matrix norm conversion error: for a CG-normalized `L -> L'` action,
> the norm summed over vector components and averaged over the source
> multiplet has weight `(2L'+1)/(2L+1)`.  The original pilot used only
> `1/(2L+1)`, underweighting target multiplets by `2L'+1`.  This especially
> changes the relative `J` and `T` channels and therefore changes both the
> fitted generator and the Hamiltonian objective.  The production output is
> now `n6_multistart_02_cg_corrected`; all numerical N=6 optimization results
> below from `n6_multistart_01` are historical diagnostics and must not be
> used as a corrected conformal-algebra result.  The same audit added an
> eigensolver retry for an exact warm vector that previously could return only
> one Ritz pair and cause a `BoundsError`.

This branch replaces a small set of hand-selected integer gap targets with a
native SO(3)lver test of a common microscopic conformal generator.  It contains
both the fixed-point audit and an experimental outer Hamiltonian optimizer.

## Representation and projection

The calculation uses only `FuzzifiED.SO3lver` composite spaces and
`CompOperator` actions.  The production profile fixes
`heavy_space_mode = "laughlin13"`, so the charge-3 segment is the exact
fermionic Laughlin-1/3 quasihole space.  No old full-Fock `L^2` projection is
performed.

The eight SU(3)-singlet rank-one candidates are the already-regressed

```text
Uf, Vf, U0, V0, Uf0, Vf0, t, mu
```

families.  A single coefficient vector is shared by the singlet and adjoint
representations.

## Algebra implemented in the pilot

For `D=(H-E0)/factor` and `Lambda=P+K`, the code applies

```text
P = (Lambda + [D,Lambda]) / 2
K = (Lambda - [D,Lambda]) / 2
```

directly between exact fixed-`L` SO(3) blocks.  Both Hamiltonian terms in the
commutator are matrix-free full-block actions.  In particular, the code does
not insert a truncated low-energy resolution of the identity between two
operators.

The first fit minimizes a generalized Rayleigh quotient containing

```text
||Lambda |vacuum>||^2 + sum_phi sum_L' [ ||K_(L->L') |phi>||^2
                                      + w_D ||([D,P]-P)|phi>||^2
                                      + w_low ||(1-Q_low)P|phi>||^2 ]
----------------------------------------------------------------
              sum_phi sum_L' ||Lambda_(L->L') |phi>||^2
```

for the selected candidate primaries `S`, `O`, `J`, and `T`.  The coefficient
normalization therefore cannot collapse to the zero generator and is
invariant under nonsingular rescalings of the candidate basis.  The cylinder
energy scale is fitted at the same time rather than fixed by declaring one
spectral gap to be an integer.  `w_D` is the explicit `dilatation_weight` in
the configuration.  The `[D,K]+K` numerator is the negative of the displayed
`[D,P]-P` numerator for `P,K` defined from the same `Lambda`, so it is not
double-counted in the fit.

Every displayed channel norm denotes the physical sum over the three vector
components, averaged over the `2L+1` magnetic substates of the source.  Since
the stored SO(3)lver matrix is CG-normalized, a reduced `L -> L'` matrix norm is
multiplied by `(2L'+1)/(2L+1)`.  The scalar `[Kz,Pz]` normalization remains a
single Cartesian component and consequently does not carry that factor of
three.

`Q_low` is the projector onto a configurable number of low-energy states in
the complete target SO(3) block.  It does not select a named descendant, but
prevents the fit from satisfying the dilation algebra by moving all generator
weight to a high-energy transition and choosing a correspondingly large
cylinder scale.

After finding the scale-invariant generator direction, the overall generator
normalization is fixed algebraically from the scalar-primary expectation value
`[Kz,Pz]=2D`.  For an `L=0,m=0` source, the reduced-matrix-element convention
actually used by the FuzzifiED builders is Clebsch--Gordan normalized: the
`L=0,m=0 -> L=1,m=0` coefficient is one.  There is therefore no extra `1/3`
in this normalization.  This was checked directly against the reverse reduced
operator, which carries the required `-sqrt(3)` adjoint factor.  The
independently normalized residuals for the scalar training primaries are
written to the primary table.  A scalar reserved as a holdout is not used for
this normalization.

The audit then reports, separately for every primary and allowed target
angular momentum:

- `||K|phi>||^2 / ||Lambda|phi>||^2`;
- the normalized `[D,P]-P` residual;
- the normalized `[D,K]+K` residual;
- `K^dagger K` eigenvalues and eigenvectors in each retained low-energy source
  subspace, allowing a finite-size primary to be a mixture of energy
  eigenstates;
- all nine Cartesian components of
  `[K_i,P_j]=2 delta_ij D-2i M_ij` on every magnetic substate of `S,O,J,T`;
- the three independent components of `[P_i,P_j]=0` and `[K_i,K_j]=0`.

Magnetic substates are reconstructed from the native reduced matrices using
the same Clebsch--Gordan convention as the FuzzifiED source.  Both orders of
each two-generator product propagate through every allowed full projected
SO(3) block.  In particular, the `T` audit includes singlet `L=4` as a final
block, but does not diagonalize it.  The Lorentz/rotation subalgebra and the
action of rotations on rank-one tensors are exact by construction; the
nontrivial finite-size tests are therefore the dilation and three commutators
listed above.

The mixed residual is reported as a squared-norm ratio
`||[K,P]-RHS||^2/||RHS||^2`, summed over Cartesian components and averaged over
the source multiplet.  Since the `P-P` and `K-K` targets vanish, their squared
norms use the same mixed-algebra RHS norm as a common scale.  These complete
commutators are currently strict audits, not terms in the generalized
eigenproblem that fits the generator direction.

## Run

```bash
julia --threads=4 --project=. scripts/so3lver_conformal_algebra.jl \
  --config=config/so3lver/n6_projected_conformal_algebra.toml
```

The strict scalar-training/spinning-holdout audit is

```bash
julia --threads=4 --project=. scripts/so3lver_conformal_algebra.jl \
  --config=config/so3lver/n6_projected_conformal_algebra_scalar_fit.toml
```

The output directory contains:

- `algebra_summary.csv`;
- `algebra_generator_coefficients.csv`;
- `algebra_primary_residuals.csv`;
- `algebra_channel_residuals.csv`;
- `algebra_mixed_commutator_residuals.csv` and
  `algebra_mixed_commutator_components.csv`;
- `algebra_k2_modes.csv` and `algebra_k2_eigenvectors.csv`;
- `algebra_metadata.toml` with the exact point, dimensions, fitted scale, and
  numerical-rank diagnostics.

Do not use `algebra_loss` in an outer Hamiltonian search until several fixed
Hamiltonian points have been audited and at least one primary family has been
reserved as a holdout.  The fit/holdout labels are explicit in the TOML and in
every result row.

## Hamiltonian optimization

`scripts/so3lver_conformal_optimize.jl` performs an actual outer search over
Hamiltonian couplings.  It does not choose `muc` afterward with the old
integer-spectrum locator.  In the supplied projected profile the four direct
optimization coordinates are

```text
Uf, Uf0, Vf0, muc
```

with `U0=9Uf`, `Vf=0`, `V0=0`, and `t=0.5`.  The exact
parameter list and bounds remain configuration data.  The reusable conformal
problem builds projected SO(3) spaces and generator operators once; every
candidate then retunes the Hamiltonian, resolves the required low-energy
states, refits the common generator, and evaluates the complete algebra.

The current training set is `S,O,J`; `T` is excluded from both the generator
fit and Hamiltonian objective.  Each training primary contributes

```text
K-primary, dilation, [K,P], [P,P], [K,K]
```

and the vacuum contributes one additional constraint.  The objective is the
weighted mean of these 16 squared-norm ratios plus `worst_weight` times their
largest value.  The state at every point must retain sufficient same-rank
overlap with the anchor primary, the generator normalization must be valid,
and the fitted cylinder factor may not lie on its configured boundary.

Run the N=6 search with

```bash
mkdir -p slurm-logs
sbatch slurm/so3lver_n6_conformal_optimize.sbatch
```

The production Slurm entry launches six independent one-CPU Julia workers.
Each worker receives one deliberately separated four-dimensional
Nelder--Mead start and four independently seeded Latin-hypercube points, for
24 exploratory LHS points in total.  Each worker may trigger at most one
automatic expansion into the wider hard bounds, and `muc` is one of the four
simultaneous simplex coordinates.  Workers have separate projected SO(3)
workspaces, CSV files, and source-signed checkpoints, so resubmitting the same
Slurm file resumes their completed Hamiltonian points without concurrent
writes.  After all six finish, a separate aggregator checks parameter and
objective consensus across their endpoints.  The corrected profile
deliberately writes a new directory,
`output/so3lver/conformal_optimization/n6_multistart_02_cg_corrected`, so it
cannot resume or overwrite the invalidly weighted `n6_multistart_01` trace.

The Slurm entry requests six tasks on one node and does not set `--mem` or
`--time`.  On this cluster the CPU request already determines the memory
allocation (currently about 7.8 GiB per CPU), so six CPUs automatically supply
ample memory for six roughly 2 GiB projected workspaces.  Each worker uses one
Julia thread; the benchmark showed no wall-time gain from adding threads
*inside one solve*, while independent multistart workers do reduce wall time.

Selection and ranking use only the `S/O/J` conformal-algebra objective.  The
final candidate is accepted only after a cold eigensolver recheck, multistart
basin agreement, coordinate-neighbor tests, the state-identity/factor gates,
and a `T` holdout veto.  The old five- and stable-six energy scores are written
to `best_spectrum_diagnostics.csv` only after selection and do not enter any
search or acceptance decision.  The main server outputs are

- `best.toml` (point, acceptance flags, cold check, holdout, and diagnostic
  spectrum scores);
- `algebra_optimization_evaluations.csv` (resumable point trace);
- per-worker `multistart_convergence.csv` and `robustness_neighbors.csv`;
- root-level `parallel_convergence.csv` and `parallel_summary.toml`;
- `top_candidates.csv` and the best-point constraint/component tables.

For a direct interactive run, use

```bash
julia --threads=1 --project=. scripts/so3lver_conformal_optimize.jl \
  --config=config/so3lver/n6_projected_conformal_optimization.toml
```

The small end-to-end regression profile is

```bash
julia --threads=4 --project=. scripts/so3lver_conformal_optimize.jl \
  --config=config/so3lver/n4_projected_conformal_optimization_smoke.toml \
  --allow-test-size=true
```

After the CG-norm correction, an N=4 smoke search reduced the training
objective from `0.421166` to `0.325437` and independently reduced the `T`
holdout mean from `0.354888` to `0.309333`; it moved `muc` from `0.228811` to
`0.255884`.  It passed the cold recheck, local-neighbor, multistart, boundary,
and holdout gates.  Before the correction, a deliberately tiny N=6 pilot (two
Latin-hypercube points and two local iterations) found

```text
Uf=0.457406, U0=4.116650, Uf0=1.707473,
Vf0=0.418754, V0=0, muc=0.276065
```

and reduced the then-misweighted training objective from `0.706920` to
`0.338621`, while the
strict `T` holdout mean fell from `0.312292` to `0.261337`.  This is evidence
only that the workflow was numerically viable, not a retained physical point:
besides the later normalization correction, the optimizer was intentionally
unconverged and its minimum anchor-state overlap was `0.69945`, just below the
production profile's `0.70` gate.  The corrected production
multistart/robustness workflow above is the test of whether a genuinely new
algebraic basin exists.

## N=6 projected audit

The first local comparison used identical algebra settings at two projected
Hamiltonian points:

| point | algebra loss | fitted factor | max primary `K^2/Lambda^2` | max `P` dilation residual | scalar `[Kz,Pz]` RMS |
|---|---:|---:|---:|---:|---:|
| latest grid point | 0.343526 | 0.086230 | 0.148090 | 0.525567 | 0.181063 |
| historical original seed | 0.345332 | 0.088738 | 0.166962 | 0.534013 | 0.157488 |

The spectrum-tuned latest point is only about 0.5% better in the combined
scale-invariant loss.  It improves the worst `K` and dilation residuals, mostly
through `T`, but is worse in the independently normalized scalar mixed
commutator.  Thus the algebra audit does not merely reproduce the old spectrum
ranking.

The complete spin-resolved commutators give a sharper and partly reversed
comparison.  The entries below are squared-norm ratios; smaller is better.

| point / fit | primary | mixed `[K,P]` | `[P,P]` | `[K,K]` |
|---|---|---:|---:|---:|
| latest grid, all-fit | `S` | 1.575518 | 0.790071 | 0.105246 |
|  | `O` | 0.511485 | 0.218866 | 0.038025 |
|  | `J` | 0.522079 | 0.195250 | 0.054469 |
|  | `T` | 0.773485 | 0.056792 | 0.015652 |
| original seed, all-fit | `S` | 1.563347 | 0.737772 | 0.097727 |
|  | `O` | 0.476355 | 0.198575 | 0.030071 |
|  | `J` | 0.493744 | 0.190626 | 0.054782 |
|  | `T` | 0.796458 | 0.053145 | 0.014374 |

Thus the latest spectrum point is better only in the `T` mixed commutator; the
original seed is better in `S`, `O`, and `J`, and usually also in translation
commutativity.  The complete-commutator audit therefore reverses the small
preference for the latest grid point in the earlier generalized fit loss.

At the latest point, fitting only `S,O` gives `factor=0.089870`.  The strict
holdouts are

| holdout | `K^2/Lambda^2` | `P` dilation residual | low-energy leakage |
|---|---:|---:|---:|
| `J` | 0.094951 | 0.423245 | 0.054072 |
| `T` | 0.169726 | 0.596291 | 0.079266 |

For this scalar-only fit the mixed residuals on `(S,O,J,T)` are respectively
`(1.548816, 0.580929, 0.515671, 0.750487)`.  Making the two scalar diagonal
expectation values nearly exact therefore does not make their complete tensor
commutator clean: the forbidden off-diagonal components remain large.

This is not yet a clean conformal point: in particular the `T` holdout remains
the dominant failure under the original primary/dilation metrics, while `S`
is the dominant failure of the full tensor algebra.  These numbers justify an
algebra-driven parameter scan, but not yet a production optimization using a
single scalar objective.  The next step is to choose fit/holdout weights for
the complete commutators and test a small Hamiltonian scan before allowing a
continuous outer optimizer to exploit them.
