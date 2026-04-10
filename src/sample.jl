"""

    _sample_mean!(d, dist, n)

Update the vector with average effects

# Arguments
 - `d`:     The databucket created with the create_databucket function
 - `dist`:  The distibution that the average effects should be sampled from
 - `n`:     The number of observations
             
# Returns
 - `Nothing`: 

"""
function _sample_mean!(d, dist, n)

    # Create aliases for temporary vectors
    random_deviations = d.tmpshared.GFloat1
    average_deviations = d.tmpshared.GFloat2
    total_deviations = d.tmpshared.GFloat3
    residuals = d.in.R

    # Sample deviations
    rand!(dist, random_deviations)

    # Calculate deviations in existing residuals
    mean!.(eachrow(average_deviations), residuals)

    # Combined observed and sampled deviations
    total_deviations .= average_deviations + sqrt.(d.out.Rvar ./ n) .* random_deviations

    # Update the average effects
    d.out.mu .+= total_deviations

    # Update the residuals
    for i in axes(total_deviations, 1)
        residuals[i] .-= total_deviations[i]
    end

    # Return
    return nothing
end

"""

   _sample_residualvariance!(d, idist)

Update the vector with residual variances

# Arguments
 - `d`:                     The databucket created with the create_databucket function
 - `idist`:                 The inverse of the distibution that the residual variances should be sampled from

 # Internal objects from the d object
 - `d.tmpshared.GFloat1`:   A vector of Float64 values
 - `d.out.Rvarscale`:       A vector of Float64 values with the prior for the residual variance
 - `d.in.R`:                A vector of vectors of Float64 values with residuals

# Returns
 - `Nothing`: 

"""
function _sample_residualvariance!(d, idist)

    # Create aliases for temporary vectors
    random_deviations = d.tmpshared.GFloat1

    # Sample deviations
    rand!(idist, random_deviations)

    # Calculate new residual variances
    d.out.Rvar .= (d.out.Rvarscale + dot.(d.in.R, d.in.R)) ./ random_deviations
    
    # Return
    return nothing
end

"""

  _sample_effectvars!(variances, effects, prior_df, n, prior_offset, random_deviations)

Update the vector with variances of gene-gene or marker-gene effects

# Arguments
 - `variances`:             A vector of Float64 values with variances for the focal effect
 - `effects`:               A vector of vectors with Float64 values for the focal effect 
 - `prior_df`:              A vector of Float64 values of the degrees of freedom for the prior variances
 - `n`:                     A vector of Float64 values 
 - `prior_offset`:          A vector of Float64 values of prior degrees of freedom times variance
 - `random_deviations`:     A vector of Float64 values for sampling random values into

# Returns
 - `Nothing`: 

"""
function _sample_effectvars!(variances, effects, prior_df, n, prior_offset, random_deviations)

    # Sample from the chi distribution into the placeholder vector
    rand!.(Chisq.(prior_df + n), eachrow(random_deviations))

    # Update the effect variances
    variances .= (prior_offset + dot.(effects, effects)) ./ random_deviations
    
    # Return
    return nothing
end

"""

  _sample_effectvars!(variances, effects, prior_df, n, prior_offset, random_deviations)

Update the vector with inclusion probabilities of gene-gene or marker-gene effects.
Both pi, n, and a are updated.

# Arguments
 - `pi`:                    A vector of Float64 values with inclusion probabilities for the focal effect
 - `indicators`:            A vector of vectors with Bool values for whether the focal effect is not zero
 - `effects`:               A vector of vectors with Float64 values for the focal effect 
 - `n`:                     A vector of Float64 values for the number of times a gene has focal effects that are not zero
 - `a`:                     A vector of Float64 values for the number of times a gene has focal effects that are zero

# Returns
 - `Nothing`: 

"""
function _sample_pi!(pi, indicators, effects, n, a)

    # Calculate the number of effects per gene
    sum!.(eachrow(n), indicators)

    # Safeguard against n=0? Not sure
    n .+= 1

    # Calculate the number of effects that are zero
    a .= ((2 + size(effects[1], 1)) .- n) # This part of the function can be optimized because it is constant across iterations.

    # Sample new pi from the Beta distribution
    rand!.(Beta.(n, a), eachrow(pi))

    # Return
    return nothing
end


function _sample_effect!(c, q, x, r, v0s, v1s, values, w, wi, Wpi, rands, a, ssqX, Rvar, norms, d, tmp)

    # Calculate intermediate values
    values.dotprod .= dot!(tmp, x, r)
    values.logD0 .= -0.5 * (log(v0s[q]) + values.dotprod[1]^2 / v0s[q]) + log(1.0 - Wpi[c])
    values.logD1 .= -0.5 * (log.(v1s[q]) + values.dotprod[1]^2 / v1s[q]) + log(Wpi[c])
    values.prob .= 1.0 / (1.0 + exp(values.logD0[1] - values.logD1[1]))


    wi[q] = rands[q] < values.prob[1]
    if wi[q] == 1
        w[q] = values.dotprod[1] / (ssqX[q] + a[c]) + sqrt(Rvar[c] / (ssqX[q] + a[c])) * norms[q]
        r .-= w[q] * x
    else
        w[q] = 0.0
    end
    return nothing
end


function _sample_effectparameters!(d, i, Z, ssqs, Normdist)

    # Aliases
    X = d.in[Z]
    ssqX = ssqs[Z]
    W = d.out[i]

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
    nsplit = min(nthreads(), size(W, 1))
    tchunks = split_vector(axes(W, 1), nsplit)

    # Work per thread
    @floop for s in 1:nsplit
        v0s = (Z == :M) ? d.tmpthread[s].QFloat1 : d.tmpthread[s].GFloat1
        v1s = (Z == :M) ? d.tmpthread[s].QFloat2 : d.tmpthread[s].GFloat2
        rands = (Z == :M) ? d.tmpthread[s].QFloat3 : d.tmpthread[s].GFloat3
        norms = (Z == :M) ? d.tmpthread[s].QFloat4 : d.tmpthread[s].GFloat4
        tmp = d.tmpthread[s].IFloat1 
        values = d.tmpthread[s].values

        # Work per gene
        for c in tchunks[s]
            v0s .= ssqX * Rvar[c]
            v1s .= ((ssqX .^ 2) * Wvar[c]) .+ v0s
            r = R[c]
            rand!(rands)
            rand!(Normdist, norms)
            w = W[c]
            wi = Wi[c]

            # Work per effect
            for q in candidates[c]
                @views x = X[:, q]
                r .+= w[q] .* x
                _sample_effect!(c, q, x, r, v0s, v1s, values, w, wi, Wpi, rands, a, ssqX, Rvar, norms, d, tmp)
            end
        end
    end

    # Sample inclusion probabilities
    _sample_pi!(Wpi, Wi, W, n, a)

    # Sample variances
	### This order is very important here, as "n" is updated in sample_pi which comes before... Could they be independent???
    _sample_effectvars!(Wvar, W, Wvardf, n, Wvarscale, a)

    # Return nothing
    return nothing
end
