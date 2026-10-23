# GRN.jl — Agent Guidelines

## Project
Julia package for Bayesian gene regulatory network inference via Gibbs sampling.
MCMC over gene-gene (G) and marker-gene (Q) effects with selectable priors.

## Architecture
- `src/GRN.jl` — module definition, imports, exports
- `src/priors.jl` — prior selectors (`SpikeSlab`, `Ridge`, `Laplace`) and state types
- `src/initialize.jl` — data bucket creation, prior insertion
- `src/sample.jl` — Gibbs samplers (allocation-free hot path)
- `src/estimate.jl` — `grn()` entry point, `_outputs()` collector
- `src/io.jl` — file output (gzip compression, IOBuffer batching, @printf formatting)
- `src/postanalysis.jl` — `accuracy`, `posterior_mean`, `posterior_inclusion_probability`, `derive_indicator`
- `src/check.jl`, `src/utility.jl` — input validation, helpers

## Key decisions (do not regress)
- **No indicator files written.** PIP is derived exactly from effects files
  via `posterior_inclusion_probability()` (counting non-zero values across rows).
- **Output is gzip-compressed by default** (`compress=true` is the default in `grn()` and `makeout()`, level 1).
  Set `compress=false` for plain `.tsv` files.
- **IO writes use IOBuffer + @printf** (`%.10g` for Float64). No per-element `print` to streams.
- **Gibbs hot path is allocation-free.** Scalar loops, in-place updates, no temporary vectors.
- **Adding a new prior**: define selector + state subtype in `priors.jl`, add
  `_initstate`, `_sample_effects_gene!`, `_sample_hyperparameters_gene!`, `_outputs`, `_varscale!`.

## Scale target
17K genes × 17K markers. Output files can be ~5 GB per iteration (uncompressed text).
Memory and write speed matter. Profile before optimizing.

## Style
- 4-space indentation, no tabs
- Type annotations on function signatures (`::Float64`, `::Nothing`, etc.)
- `@inbounds` in hot loops
- Docstrings with `"""` for all public and internal functions
- NamedTuples for structured data (no mutable structs for state)

## Testing
- Run with `Pkg.test()` (from package directory)
- Tests in `test/runtests.jl` cover: input checks, prior insertion, scalar samplers,
  output formats for all priors, compressed output, PIP derivation, allocations
- Add tests for any new output file or post-analysis function

## Dependencies
- Standard libs: `Printf`, `LinearAlgebra`, `Random`, `Statistics`, `DelimitedFiles`, `GZip`
- No `Distributions` — all random draws use internal scalar samplers
- When adding a dependency, update `Project.toml` `[deps]` and `[compat]`,
  then run `Pkg.resolve()`

## Changelog
- **Every change must be recorded in `CHANGELOG.md`.** Update it in the same
  commit as the code change — never in a separate commit.
- Use the existing `[Unreleased]` section at the top. Do not create a new
  version header unless releasing.
- Categorize entries under the existing subsections:
  - `### Added` — new features, functions, files, priors, options
  - `### Changed` — behavior changes, refactors, signature changes, defaults
  - `### Removed` — deleted code, dependencies, output files
  - `### Performance` — speed or memory improvements (include before/after numbers)
- Write entries as bullet points with a short description. Reference the file
  or function affected, e.g. `` `grn()` now accepts `compress=true` (`src/estimate.jl`) ``.
- If a change breaks backward compatibility, say so explicitly at the top of
  the entry, e.g. "**Breaking:** `grn()` no longer writes `G_indicator.tsv`."
