# Full-algebra implementation validation, 2026-09-27

Environment: Julia 1.11.6, FuzzifiED 2.0.1, one Julia/BLAS thread per process.
All numerical execution here is small-size (N2/N4); no N6/N7 ED was run.

## Operator and algebra regressions

- `test/so3lver_runtests.jl`: 216/216 pass, including spectra versus the
  conventional construction, exact Laughlin projection and CG norm weighting.
- `test/so3lver_conformal_fit_runtests.jl`: 179/179 pass:
  9 exact free-scalar algebra checks including spinning descendants;
  6 continuation/level-crossing checks; 123 projected density-basis checks;
  41 inner/outer objective, holdout and noncollapse checks.
- The missing light-Laplacian/heavy density moment has relative projection
  residual greater than 0.70 in the old eight-candidate basis, below 1e-10
  in the new 15-candidate basis. Candidate Hermiticity is tested independently
  across both SU(3) representations and four angular-momentum transitions.
- Compact full-commutator formulas agree with direct operator application,
  including perturbed coefficients and scale, for all five training sources.
  A separate diagnostic obtained differences around 1e-13; assertions use
  1e-8 tolerance. Every production evaluation repeats the final-objective
  comparison at 1e-7 tolerance.
- At the N4 anchor with six-state leakage projectors, the same full objective
  drops from 0.7987703358711105 at the proxy initializer to
  0.3133305019654563 after joint refinement; the inner optimizer converges.
  These are two generators at the **same Hamiltonian**, not a new CFT point.

## Workflow tests

Two real independent N4 workers ran simultaneously with distinct seeds,
one exploratory sample, one local start, and one outer iteration each.
The N4 smoke profile intentionally has loose acceptance/consensus thresholds;
its acceptance flag is an execution-path test, not a scientific claim.

| quantity | worker 1 | worker 2 |
|---|---:|---:|
| anchor training objective | 0.3185416712099201 | 0.3185416712099201 |
| selected training objective | 0.28616446253765304 | 0.31028804268944526 |
| selected T holdout | 0.25360372197517356 | 0.33030458592202033 |
| cold-recheck difference | 0 | 0 |

The smoke profile uses four-state leakage projectors, hence its anchor
objective differs from the six-state operator regression above.

Both workers finish coefficient/residual export, eight coordinate-neighbor
checks, and final metadata. The real two-worker aggregator runs successfully
and checks matching source/configuration hashes and inner convergence.
Resuming worker 1 reuses its eight search points and excludes nine prior
final-audit points from search selection; only the final audits are repeated.
An incompatible older output signature is rejected before building spaces.

The fixed-point pilot was also exercised at N4 with a second point
`(Uf,Uf0,Vf0,mu)=(0.49,1.70,0.36,0.235)`. All three inner fits converge;
cold repeat difference is zero. `pilot.toml`, residuals, coefficients and trace
are emitted. N6 configuration validation and both Slurm syntax checks pass.

## Remaining limits

N6 convergence, cost, and any improvement over previous Hamiltonian points
require the server pilot. An exact free-scalar module is a kinematic positive
control, not an interacting fuzzy-Ising benchmark. The finite generator
ansatz, local coefficient optimization, finite-size corrections, rank-based
state hypotheses, and full-energy P/K UV sensitivity remain limitations.
Do not interpret a lower loss alone as evidence of the desired CFT.

See `SO3LVER_CONFORMAL_V3_SERVER.md` for update, one-CPU pilot, six-worker
search, and download instructions.
