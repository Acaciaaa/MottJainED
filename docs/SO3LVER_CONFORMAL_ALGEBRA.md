# Projected conformal-algebra pilot

This branch replaces a small set of hand-selected integer gap targets with a
native SO(3)lver test of a common microscopic conformal generator.  It is an
audit/prototyping path, not yet the production outer-parameter optimizer.

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

## Algebra used in the first pilot

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

`Q_low` is the projector onto a configurable number of low-energy states in
the complete target SO(3) block.  It does not select a named descendant, but
prevents the fit from satisfying the dilation algebra by moving all generator
weight to a high-energy transition and choosing a correspondingly large
cylinder scale.

After finding the scale-invariant generator direction, the overall generator
normalization is fixed algebraically from the scalar-primary expectation value
`[Kz,Pz]=2D`.  For an `L=0,m=0` source, the reduced-matrix-element convention
used by SO3lver supplies the `1/3` Wigner-Eckart factor.  The independently
normalized residuals for the scalar training primaries are written to the
primary table.  A scalar reserved as a holdout is not used for this
normalization.

The audit then reports, separately for every primary and allowed target
angular momentum:

- `||K|phi>||^2 / ||Lambda|phi>||^2`;
- the normalized `[D,P]-P` residual;
- the normalized `[D,K]+K` residual;
- `K^dagger K` eigenvalues and eigenvectors in each retained low-energy source
  subspace, allowing a finite-size primary to be a mixture of energy
  eigenstates.

Rotational covariance is exact by construction because the operators are
native rank-one SO(3) tensors.  A coupled-tensor implementation of the mixed
commutator `[K_i,P_j]=2 delta_ij D-2i M_ij`, followed by `[P_i,P_j]` and
`[K_i,K_j]`, remains a separate second step; it must retain the SO(3) reduced
matrix-element normalization and all full intermediate blocks.

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
- `algebra_k2_modes.csv` and `algebra_k2_eigenvectors.csv`;
- `algebra_metadata.toml` with the exact point, dimensions, fitted scale, and
  numerical-rank diagnostics.

Do not use `algebra_loss` in an outer Hamiltonian search until several fixed
Hamiltonian points have been audited and at least one primary family has been
reserved as a holdout.  The fit/holdout labels are explicit in the TOML and in
every result row.

## First N=6 projected audit

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

At the latest point, fitting only `S,O` gives `factor=0.089870`.  The strict
holdouts are

| holdout | `K^2/Lambda^2` | `P` dilation residual | low-energy leakage |
|---|---:|---:|---:|
| `J` | 0.094951 | 0.423245 | 0.054072 |
| `T` | 0.169726 | 0.596291 | 0.079266 |

This is not yet a clean conformal point: in particular the `T` holdout remains
the dominant failure.  These numbers justify an algebra-driven parameter scan,
but not yet a production optimization using a single scalar objective.  The
next implementation step should add the spinning-state components of
`[K_i,P_j]=2 delta_ij D-2iM_ij` and compare several more fixed Hamiltonian
points before fixing the objective weights.
