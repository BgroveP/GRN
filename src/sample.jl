

function _sample_residualvariance!(d, idist)
    rand!(idist, d.tmpshared.GFloat1)
    d.out.Rvar .= (d.out.Rvarscale + dot.(eachcol(d.in.R), eachcol(d.in.R))) ./ d.tmpshared.GFloat1
    return nothing
end


function _sample_mean!(d, dist, n)

    # Create aliases
    f1 = d.tmpshared.GFloat1
    f2 = d.tmpshared.GFloat2
    f3 = d.tmpshared.GFloat3
    r = d.in.R

    # Work
    rand!(dist, f1)
    mean!.(eachrow(f2), eachcol(r))
    f3 .= f2 + sqrt.(d.out.Rvar ./ n) .* f1
    d.out.mu .+= f3
    r .-= f3'
    return nothing
end

"""
W is the coefficient matrix (you called in Lambda)

"""
######### PLEASE FIND A DIFFERENT NAME. YOU ALSO HAVE _sample_effect :) ##########
function _sample_effects!(d, i, Z, ssqs, Normdist; omitdiag=false)

    # Aliases
    X = d.in[Z]
    ssqX = ssqs[Z]
    W = d.out[i]
	######  WHEN DO YOU UPDATE n ######
    n = d.tmpshared.GFloat1 .= 0
    Rvar = d.out.Rvar
    a = d.tmpshared.GFloat2 .= Rvar ./ d.out[Symbol.(String(i) .* "var")]
    R = d.in.R
    Wvar = d.out[Symbol(string(i) * "var")]
    Wvardf = d.out[Symbol(string(i) * "vardf")]
    Wvarscale = d.out[Symbol(string(i) * "varscale")]
    Wi = d.out[Symbol(string(i) * "i")]
    Wpi = d.out[Symbol(string(i) * "pi")]
    candidates = d.out.candidates[i]

    # Prepare for multithreading
    nsplit = min(nthreads(), size(W, 2))
    tchunks = split_vector(axes(W, 2), nsplit)

    # Work per thread
    @threads for s in 1:nsplit
        v0s = (Z == :M) ? d.tmpthread[s].QFloat1 : d.tmpthread[s].GFloat1
        v1s = (Z == :M) ? d.tmpthread[s].QFloat2 : d.tmpthread[s].GFloat2
        rands = (Z == :M) ? d.tmpthread[s].QFloat3 : d.tmpthread[s].GFloat3
        norms = (Z == :M) ? d.tmpthread[s].QFloat4 : d.tmpthread[s].GFloat4
        values = d.tmpthread[s].values
        rng = d.tmpthread[s].rng

        # Work per gene
        for c in tchunks[s]
            v0s .= ssqX * Rvar[c]
            v1s .= ((ssqX .^ 2) * Wvar[c]) .+ v0s
            @views r = R[:, c]
            rand!(rng, rands)
            rand!(rng, Normdist, norms)
            @views w = W[:, c]
            @views wi = Wi[:, c]

            # Work per effect
            for q in candidates[c]
                @views x = X[:, q]
                @lock d.tmpshared.wlock r .+= w[q] .* x
                if omitdiag
                    if ~isequal(c,q)
                        _sample_effect!(c, q, x, r, v0s, v1s, values, w, wi, Wpi, rands, a, ssqX, Rvar, norms, d)
                    end
                else
                    _sample_effect!(c, q, x, r, v0s, v1s, values, w, wi, Wpi, rands, a, ssqX, Rvar, norms, d)
                end
            end
        end
    end

    # Sample inclusion probabilities
    _sample_pi(Wi, Wpi, W, n, a)

    # Sample variances
    _sample_effectvars(Wvar, W, Wvardf, n, Wvarscale, a)

    # Return nothing
    return nothing
end

function _sample_effectvars(Wvar, W, Wvardf, n, Wvarscale, a)
    rand!.(Chisq.(Wvardf + n), eachrow(a))
    Wvar .= (Wvarscale + dot.(eachcol(W), eachcol(W))) ./ a
    return nothing
end

function _sample_pi(Wi, Wpi, W, n, a)
    sum!.(eachrow(n), eachcol(Wi))
    n .+= 1
    a .= ((2 + size(W, 1)) .- n)
    rand!.(Beta.(n, a), eachrow(Wpi))
    return nothing
end

function _sample_effect!(c, q, x, r, v0s, v1s, values, w, wi, Wpi, rands, a, ssqX, Rvar, norms, d)

    # Calculate intermediate values
    values.dotprod .= dot(x, r)
    values.logD0 .= -0.5 * (log(v0s[q]) + values.dotprod[1]^2 / v0s[q]) + log(1.0 - Wpi[c])
    values.logD1 .= -0.5 * (log.(v1s[q]) + values.dotprod[1]^2 / v1s[q]) + log(Wpi[c])
    values.prob .= 1.0 / (1.0 + exp(values.logD0[1] - values.logD1[1]))


    wi[q] = rands[q] < values.prob[1]
    if wi[q]
        w[q] = values.dotprod[1] / (ssqX[q] + a[c]) + sqrt(Rvar[c] / (ssqX[q] + a[c])) * norms[q]
        @lock d.tmpshared.wlock r .-= w[q] * x
    else
        w[q] = 0.0
    end
    return nothing
end
