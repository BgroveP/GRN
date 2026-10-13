# Effect priors for gene-gene (G) and marker-gene (Q) effects
#
# A prior consists of
#  1. A selector type that the user constructs (e.g. `SpikeSlab()`),
#     carrying default settings for the hyperprior
#  2. A state type holding the per-gene Gibbs state, created with `_initstate`
#  3. Methods for `_sample_effects_gene!`, `_sample_hyperparameters!`,
#     `_outputs`, and `_varscale!`
#
# Adding a new prior: define a selector and a state subtype of
# `AbstractEffectPriorState`, plus the four methods above.

abstract type AbstractEffectPrior end
abstract type AbstractEffectPriorState end

# ---------------------------------------------------------------- Selectors

"""

    SpikeSlab(; includeprob, var, vardf)

Spike-and-slab prior for gene-gene or marker-gene effects. Each effect is
either zero (spike) with probability `1 - pi` or drawn from a normal
distribution with variance `var` (slab). `pi` and `var` are sampled per gene
with `vardf` degrees of freedom for the scaled inverse chi-squared prior on
`var`. Construct as a value in the `priors` NamedTuple, e.g.
`priors = (Gprior = SpikeSlab(),)`.

"""
struct SpikeSlab <: AbstractEffectPrior
    includeprob::Float64
    var::Float64
    vardf::Float64
end
SpikeSlab(; includeprob::Float64=0.5, var::Float64=1.0, vardf::Float64=4.0) =
    SpikeSlab(includeprob, var, vardf)

"""

    Ridge(; var, vardf)

Gaussian prior for all effects of a gene (ridge-like shrinkage, no spike).
The prior variance `var` is sampled per gene with `vardf` degrees of freedom
for the scaled inverse chi-squared prior. Construct as a value in the
`priors` NamedTuple, e.g. `priors = (Gprior = Ridge(),)`.

"""
struct Ridge <: AbstractEffectPrior
    var::Float64
    vardf::Float64
end
Ridge(; var::Float64=1.0, vardf::Float64=4.0) = Ridge(var, vardf)

"""

    Laplace(; lambda2, a, b)

Bayesian lasso prior (Laplace) via a scale mixture of normals: each effect
has its own local variance `tau2` with an exponential prior with rate
`lambda2 / 2`, and `lambda2` per gene has a `Gamma(a, b)` (shape, rate)
hyperprior. Construct as a value in the `priors` NamedTuple, e.g.
`priors = (Gprior = Laplace(),)`.

"""
struct Laplace <: AbstractEffectPrior
    lambda2::Float64
    a::Float64
    b::Float64
end
Laplace(; lambda2::Float64=1.0, a::Float64=1.0, b::Float64=1e-3) =
    Laplace(lambda2, a, b)

# ---------------------------------------------------------------- States

struct SpikeSlabState <: AbstractEffectPriorState
    effects::Vector{Vector{Float64}}
    indicator::Vector{Vector{Int8}}
    pi::Vector{Float64}
    var::Vector{Float64}
    vardf::Vector{Float64}
    varscale::Vector{Float64}
    candidates::Vector{Vector{Int}}
end

struct RidgeState <: AbstractEffectPriorState
    effects::Vector{Vector{Float64}}
    var::Vector{Float64}
    vardf::Vector{Float64}
    varscale::Vector{Float64}
    candidates::Vector{Vector{Int}}
end

struct LaplaceState <: AbstractEffectPriorState
    effects::Vector{Vector{Float64}}
    tau2::Vector{Vector{Float64}}
    lambda2::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    candidates::Vector{Vector{Int}}
end

# ---------------------------------------------------------------- Init

function _candidates(genes::Int, predictors::Int; diag::Bool=true)::Vector{Vector{Int}}
    o = [collect(1:predictors) for _ in 1:genes]
    if !diag
        for c in 1:genes
            deleteat!(o[c], c)
        end
    end
    return o
end

function _initstate(sel::SpikeSlab, genes::Int, predictors::Int; diag::Bool=true)
    return SpikeSlabState(
        [zeros(predictors) for _ in 1:genes],
        [zeros(Int8, predictors) for _ in 1:genes],
        fill(sel.includeprob, genes),
        fill(sel.var, genes),
        fill(sel.vardf, genes),
        fill(sel.var * (sel.vardf - 2.0), genes),
        _candidates(genes, predictors; diag=diag))
end

function _initstate(sel::Ridge, genes::Int, predictors::Int; diag::Bool=true)
    return RidgeState(
        [zeros(predictors) for _ in 1:genes],
        fill(sel.var, genes),
        fill(sel.vardf, genes),
        fill(sel.var * (sel.vardf - 2.0), genes),
        _candidates(genes, predictors; diag=diag))
end

function _initstate(sel::Laplace, genes::Int, predictors::Int; diag::Bool=true)
    return LaplaceState(
        [zeros(predictors) for _ in 1:genes],
        [fill(1.0 / sel.lambda2, predictors) for _ in 1:genes],
        fill(sel.lambda2, genes),
        fill(sel.a, genes),
        fill(sel.b, genes),
        _candidates(genes, predictors; diag=diag))
end

# ---------------------------------------------------------------- Outputs

_outputs(st::SpikeSlabState) = (
    effects=st.effects,
    indicator=st.indicator,
    pi=st.pi,
    variance=st.var)

_outputs(st::RidgeState) = (
    effects=st.effects,
    variance=st.var)

_outputs(st::LaplaceState) = (
    effects=st.effects,
    tau2=st.tau2,
    lambda2=st.lambda2)

# Prior variance scale for the scaled inverse chi-squared variance priors
_varscale!(st::SpikeSlabState) = st.varscale .= st.var .* (st.vardf .- 2.0)
_varscale!(st::RidgeState) = st.varscale .= st.var .* (st.vardf .- 2.0)
_varscale!(st::LaplaceState) = nothing
