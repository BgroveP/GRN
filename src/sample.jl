# Gibbs samplers for the gene regulatory network model.
#
# All samplers in this file are written to be allocation-free: they operate
# in place on the databucket and use only scalar operations and scalar
# random draws. Random draws for the chi-squared, beta, and gamma
# distributions are generated with a scalar version of the Marsaglia and
# Tsang (2000) gamma sampler ("A Simple Method for Generating Gamma
# Variables", ACM TOMS) to avoid the per-draw allocations of
# distribution objects.

# ---------------------------------------------------------------- Random draws

"""

    _rand_gamma(alpha)

Draw from a Gamma(alpha, 1) distribution with the Marsaglia and Tsang
(2000) squeeze-rejection method ("A Simple Method for Generating Gamma
Variables", ACM TOMS). Scalar operations only, no allocations.

"""
function _rand_gamma(alpha::Float64)::Float64
    if alpha < 1.0
        # Boost by one and multiply with a power of a uniform draw
        u = rand()
        u <= 0.0 && return 0.0
        return _rand_gamma(alpha + 1.0) * u^(1.0 / alpha)
    end
    d = alpha - 1.0 / 3.0
    c = 1.0 / sqrt(9.0 * d)
    while true
        x = randn()
        v = (1.0 + c * x)^3
        v <= 0.0 && continue
        u = rand()
        x2 = x * x
        if u < 1.0 - 0.0331 * x2 * x2
            return d * v
        end
        if log(u) < 0.5 * x2 + d * (1.0 - v + log(v))
            return d * v
        end
    end
end

_rand_chisq(df::Float64)::Float64 = 2.0 * _rand_gamma(df / 2.0)

function _rand_beta(a::Float64, b::Float64)::Float64
    x = _rand_gamma(a)
    y = _rand_gamma(b)
    return x / (x + y)
end

_rand_gamma(alpha::Float64, scale::Float64)::Float64 = scale * _rand_gamma(alpha)

"""

    _rand_invgauss(mu, lambda)

Draw from an InverseGaussian(mu, lambda) distribution with the Michael,
Schucany, and Haas (1976) algorithm. Scalar operations only, no allocations.

"""
function _rand_invgauss(mu::Float64, lambda::Float64)::Float64
    y = randn()
    y *= y
    x = mu + mu * mu * y / (2.0 * lambda) -
        mu / (2.0 * lambda) * sqrt(4.0 * mu * lambda * y + mu * mu * y * y)
    x = x > 0.0 ? x : mu
    return rand() <= mu / (mu + x) ? x : mu * mu / x
end

"""

    _sample_iteration!(d, ssqs, individuals, loci)

Perform one Gibbs iteration over all genes, in parallel over contiguous
chunks of genes.

"""
function _sample_iteration!(d::NamedTuple, ssqs::NamedTuple,
    individuals::Int, loci::Int)::Nothing

    # Aliases
    residual_variance::Vector{Float64} = d.out.base.residual_variance
    G = d.out.G
    Q = d.out.Q
    Xg::Vector{Vector{Float64}} = d.in.Y
    Xm::Vector{Vector{Float64}} = d.in.M
    R::Vector{Vector{Float64}} = d.in.R
    ssqg::Vector{Float64} = ssqs.Y
    ssqm::Vector{Float64} = ssqs.M
    chunks::Vector{Vector{Int}} = d.tmpshared.split_chunks

    # Work per parallel task
    @threads for s in eachindex(chunks)
        @inbounds for c in chunks[s]

            # Base: fixed effects (intercept is the coefficient of the ones column)
            _sample_fixedeffects!(d, c)
            _sample_residualvariance!(d, c)

            # Gene-gene effects
            _sample_effects_gene!(G, c, Xg, ssqg, R, residual_variance)
            _sample_hyperparameters_gene!(G, c, residual_variance)

            # Marker-gene effects
            if loci > 0
                _sample_effects_gene!(Q, c, Xm, ssqm, R, residual_variance)
                _sample_hyperparameters_gene!(Q, c, residual_variance)
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------- Fixed effects

"""

    _sample_fixedeffects!(d, c)

Update the fixed effect coefficients for gene `c`. With a flat prior,
each coefficient is sampled coordinate-wise from its conditional
posterior: Beta_k | rest ~ Normal(X_k'R / X_k'X_k, Rvar / X_k'X_k).
The residual is adjusted in place. When `X` is a column of ones
(the default), the single coefficient is the per-gene intercept.
Allocation-free.

"""
function _sample_fixedeffects!(d::NamedTuple, c::Int)::Nothing
    X::Vector{Vector{Float64}} = d.in.X
    beta::Vector{Float64} = d.out.base.fixedeffects[c]
    r::Vector{Float64} = d.in.R[c]
    rvar::Float64 = d.out.base.residual_variance[c]
    ssqX::Vector{Float64} = d.tmpshared.ssqX

    @inbounds for k in eachindex(X)
        x::Vector{Float64} = X[k]
        bk::Float64 = beta[k]
        prec::Float64 = ssqX[k]
        r .+= bk .* x
        bk = dot(x, r) / prec + sqrt(rvar / prec) * randn()
        beta[k] = bk
        r .-= bk .* x
    end

    return nothing
end

# ---------------------------------------------------------------- Residual variance

"""

    _sample_residualvariance!(d, c)

Update the residual variance of gene `c` from the scaled inverse
chi-squared posterior: Rvar | R ~ scale-chi2(Rvarscale + R'R, df + n).

"""
function _sample_residualvariance!(d::NamedTuple, c::Int)::Nothing
    r = d.in.R[c]
    d.out.base.residual_variance[c] = (d.out.base.residual_variancescale[c] + dot(r, r)) /
        _rand_chisq(d.out.base.residual_df[c] + length(r))
    return nothing
end

# ---------------------------------------------------------------- Effects

# The per-effect Gibbs update has the same structure for all priors:
#  1. Add the current effect back into the residual
#  2. Compute the dot product between predictor and residual
#  3. Draw the new effect from its conditional posterior
#  4. Remove the new effect from the residual
# Only step 3 differs between priors, and the shared steps are inlined
# in each specialization below to keep the hot loop type-stable.

# Spike-and-slab: sample the indicator from the mixture posterior, then the
# effect from the slab if included.
#
# The log-odds logD0 - logD1 simplifies algebraically: log(ssqX[q]) appears
# in both log(v0) and log(v1) and cancels, leaving one log call per effect
# instead of four.  log(1-pic), log(pic), and rvar/varc are hoisted out of
# the loop as per-gene constants.
function _sample_effects_gene!(st::SpikeSlabState, c::Int,
    X::Vector{Vector{Float64}}, ssqX::Vector{Float64},
    R::Vector{Vector{Float64}}, Rvar::Vector{Float64})::Nothing

    # Aliases
    w::Vector{Float64} = st.effects[c]
    wi::Vector{Int8} = st.indicator[c]
    r::Vector{Float64} = R[c]
    rvar::Float64 = Rvar[c]
    varc::Float64 = st.var[c]
    pic::Float64 = st.pi[c]
    candidates::Vector{Int} = st.candidates[c]

    # Per-gene constants hoisted out of the per-effect loop
    log_odds_prior::Float64 = log(1.0 - pic) - log(pic)  # log((1-pic)/pic)
    rvar_over_varc::Float64 = rvar / varc

    # Work per effect
    @inbounds for q in candidates
        x::Vector{Float64} = X[q]
        wq::Float64 = w[q]
        wq != 0.0 && (r .+= wq .* x)
        dotprod::Float64 = dot(x, r)
        dp2::Float64 = dotprod * dotprod
        ssqx::Float64 = ssqX[q]
        v0::Float64 = ssqx * rvar
        v1::Float64 = ssqx * (ssqx * varc + rvar)
        # logD0 - logD1 = -0.5 * (log(v0/v1) + dp2*(1/v0 - 1/v1)) + log((1-pic)/pic)
        # v0/v1 = rvar / (ssqx*varc + rvar)
        log_diff::Float64 = -0.5 * (log(rvar / (ssqx * varc + rvar)) +
            dp2 * (1.0 / v0 - 1.0 / v1)) + log_odds_prior
        wi[q] = rand() < 1.0 / (1.0 + exp(log_diff))
        if wi[q] == 1
            prec::Float64 = ssqx + rvar_over_varc
            wq = dotprod / prec + sqrt(rvar / prec) * randn()
            w[q] = wq
            r .-= wq .* x
        else
            w[q] = 0.0
        end
    end

    return nothing
end

# Ridge: all effects included, conjugate normal conditional posterior.
function _sample_effects_gene!(st::RidgeState, c::Int,
    X::Vector{Vector{Float64}}, ssqX::Vector{Float64},
    R::Vector{Vector{Float64}}, Rvar::Vector{Float64})::Nothing

    # Aliases
    w::Vector{Float64} = st.effects[c]
    r::Vector{Float64} = R[c]
    rvar::Float64 = Rvar[c]
    precision::Float64 = rvar / st.var[c]
    candidates::Vector{Int} = st.candidates[c]

    # Work per effect
    @inbounds for q in candidates
        x::Vector{Float64} = X[q]
        wq::Float64 = w[q]
        wq != 0.0 && (r .+= wq .* x)
        prec::Float64 = ssqX[q] + precision
        wq = dot(x, r) / prec + sqrt(rvar / prec) * randn()
        w[q] = wq
        r .-= wq .* x
    end

    return nothing
end

# Laplace: scale mixture of normals with per-effect local variance tau2.
# The conditional posterior of the effect is normal with precision
# ssqX + 1 / tau2 (prior variance Rvar * tau2).
function _sample_effects_gene!(st::LaplaceState, c::Int,
    X::Vector{Vector{Float64}}, ssqX::Vector{Float64},
    R::Vector{Vector{Float64}}, Rvar::Vector{Float64})::Nothing

    # Aliases
    w::Vector{Float64} = st.effects[c]
    tau2::Vector{Float64} = st.tau2[c]
    r::Vector{Float64} = R[c]
    rvar::Float64 = Rvar[c]
    candidates::Vector{Int} = st.candidates[c]

    # Work per effect
    @inbounds for q in candidates
        x::Vector{Float64} = X[q]
        wq::Float64 = w[q]
        wq != 0.0 && (r .+= wq .* x)
        prec::Float64 = ssqX[q] + 1.0 / tau2[q]
        wq = dot(x, r) / prec + sqrt(rvar / prec) * randn()
        w[q] = wq
        r .-= wq .* x
    end

    return nothing
end

# ---------------------------------------------------------------- Hyperparameters

# Spike-and-slab: sample the inclusion probability from the beta posterior
# and the prior effect variance from the scaled inverse chi-squared
# posterior, per gene.
function _sample_hyperparameters_gene!(st::SpikeSlabState, c::Int, ::Vector{Float64})::Nothing

    # Aliases
    w::Vector{Float64} = st.effects[c]
    wi::Vector{Int8} = st.indicator[c]
    n::Float64 = sum(wi)

    # Sample the inclusion probability and the prior effect variance
    st.pi[c] = _rand_beta(n + 1.0, length(w) - n + 1.0)
    st.var[c] = (st.varscale[c] + dot(w, w)) / _rand_chisq(st.vardf[c] + n)

    return nothing
end

# Ridge: sample the prior effect variance with all effects counted as
# included.
function _sample_hyperparameters_gene!(st::RidgeState, c::Int, ::Vector{Float64})::Nothing

    w::Vector{Float64} = st.effects[c]
    st.var[c] = (st.varscale[c] + dot(w, w)) /
        _rand_chisq(st.vardf[c] + length(st.candidates[c]))

    return nothing
end

# Laplace (Park and Casella, 2008): sample the local variances tau2 from
# an inverse Gaussian distribution and the penalty lambda2 per gene from a
# gamma distribution. The residual variances enter the conditionals, so the
# prior variance of an effect is Rvar * tau2.
function _sample_hyperparameters_gene!(st::LaplaceState, c::Int,
    Rvar::Vector{Float64})::Nothing

    # Aliases
    w::Vector{Float64} = st.effects[c]
    tau2::Vector{Float64} = st.tau2[c]
    rvar::Float64 = Rvar[c]
    lambda2::Float64 = st.lambda2[c]
    candidates::Vector{Int} = st.candidates[c]

    # Sample the local variances
    s::Float64 = 0.0
    @inbounds for q in candidates
        wq::Float64 = w[q]
        s += tau2[q]
        if abs(wq) > 1e-300
            # tau2 | w, lambda2, Rvar is 1 / InverseGaussian
            mu::Float64 = sqrt(lambda2 * rvar / (wq * wq))
            tau2[q] = 1.0 / _rand_invgauss(mu, lambda2)
        else
            # Limiting case: Gamma(1/2, rate lambda2 / 2)
            tau2[q] = _rand_gamma(0.5, 2.0 / lambda2)
        end
    end

    # Sample the penalty: lambda2 | tau2 ~ Gamma(a + p, rate = b + sum(tau2) / 2)
    st.lambda2[c] = _rand_gamma(st.a[c] + length(candidates), 1.0 / (st.b[c] + s / 2.0))

    return nothing
end
