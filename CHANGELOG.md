# Changelog

## [Unreleased] - 2026-10-13

Rewrite of the sampler core: the Gibbs hot path is now allocation-free, and
the prior for gene-gene and marker-gene effects is selectable per relation
type instead of being fixed to spike-and-slab.

### Added

- Prior selection with selector types in the `priors` NamedTuple
  (`src/priors.jl`):
  - `SpikeSlab(; includeprob, var, vardf)` - spike-and-slab prior (default,
    previous behavior)
  - `Ridge(; var, vardf)` - Gaussian prior on all effects with a
    sampled per-gene prior variance
  - `Laplace(; lambda2, a, b)` - Bayesian lasso prior with per-effect local
    variances (scale mixture of normals, Park and Casella 2008)
  - Selected with `priors = (Gprior = ..., Qprior = ...)`; a new prior is
    added as a selector and state subtype plus three methods
- Hyperparameters for the new priors can be set with the `priors`
  NamedTuple: `Gtau2`, `Glambda2`, `Ga`, `Gb` (and the `Q` equivalents)
- Convenience method `grn(Y, chain; ...)` for analyses without genotypes
- Marker-gene output files are no longer created when there are no marker
  loci
- Input checks for the `grn` function (`src/check.jl`), which is now
  included in the module: matching row counts, finite values, non-empty
  increasing chain, and valid prior selectors
- Test suite (`test/runtests.jl`): input checks, prior insertion, scalar
  samplers, output formats for all priors, network recovery on simulated
  data, and allocation assertions
- Benchmark scripts in `exclude/` (not part of the package): `bench.jl`
  (timing, allocations, and recovery on the data in `data/`) and
  `allocbench.jl` (per-function allocation counts)

### Changed

- `grn` runs: the main loop previously called undefined variables and
  `create_databucket` referenced undefined names, so the package did not
  run; all wiring is fixed and verified end to end
- Allocation-free Gibbs hot path:
  - Per-effect and per-hyperparameter updates are scalar loops in place
    on the data bucket; no temporary vectors are created
  - The whole per-gene sweep (mean, residual variance, gene-gene effects,
    marker-gene effects, and hyperparameters) is fused into a single
    threaded region per iteration, since each stage depends only on the
    residual of its own gene
  - Chi-squared, beta, and gamma draws use a scalar Marsaglia and Tsang
    (2000) gamma sampler; the Laplace local variances use the Michael,
    Schucany, and Haas (1976) inverse Gaussian sampler
  - The remaining allocation per iteration is the task creation of
    `@threads` (about 0.9 kB per task); a persistent task pool would be
    needed to remove it
  - Measured on `data/` (1305 individuals, 10 genes, 200 markers,
    12 threads): 0 bytes for the Gibbs updates, about 8 kB per iteration
    total, 0.28-0.41 ms per iteration
- The sweep order changed from staged (all means, then all residual
  variances, then all gene-gene effects, then all marker-gene effects) to
  per-gene fused. Both are valid Gibbs orderings, but results are not
  identical to previous runs even with the same seed
- `insert_priors!` maps keys to fields by name (base fields by field name,
  prior fields prefixed with the group letter) and recalculates the prior
  variance scales unless a `varscale` key is given explicitly
- Output files are created from a flat list of output objects: one
  `<group>_<parameter>.tsv` file per object, with names that depend on the
  selected prior (for example `G_tau2.tsv` and `G_lambda2.tsv` only for the
  Laplace prior)
- Residual variance defaults: the prior residual variance scale is now
  initialized consistently as `Rvar * (Rvardf - 2)`

### Removed

- `Distributions` dependency: all draws are generated with scalar samplers
  internal to the package
- Internal helpers that are no longer used: `_sample_pi!`,
  `_sample_effectvars!`, `dot!`, `indicator2indice`, and the temporary
  per-thread buffers (`GFloat*`, `QFloat*`, `IFloat*`) from the data bucket

### Performance

Comparison on the simulated data in `data/` (10000 iterations with 5000
burn-in, posterior mean over 51 kept samples, recovery measured as the
correlation with the true effects in `B.csv` and `F.csv`):

| Prior      | Time  | cor(G) | cor(Q) |
|------------|-------|--------|--------|
| SpikeSlab  | 5.9 s | 0.770  | 0.468  |
| Laplace    | 12.4 s | 0.805  | 0.518  |
| Ridge      | 10.2 s | 0.739  | 0.444  |
