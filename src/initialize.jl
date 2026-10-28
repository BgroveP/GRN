"""

    create_databucket(Y, M, amounts, Gsel, Qsel, X=ones(n, 1))

Create the data container with input data, parameters, and prior states.

# Arguments
 - `Y`:             An *n* × *m* matrix with gene expressions   (Float64)
 - `M`:             An *n* × *u* matrix with genotypes          (Float64)
 - `amounts`:       A NamedTuple with individuals, genes, and loci  (Int)
 - `Gsel`:          The prior selector for gene-gene effects    (AbstractEffectPrior)
 - `Qsel`:          The prior selector for marker-gene effects (AbstractEffectPrior)
 - `X`:             An *n* × *p* matrix with fixed effect covariates.
                    Defaults to `ones(n, 1)` (intercept only). The per-gene
                    mean is the coefficient of the intercept column.

# Returns
 - `NamedTuple`: (in = (; ...), out = (base = (; ...), G = <prior state>, Q = <prior state>))

"""
function create_databucket(Y::Matrix{Float64}, M::Matrix{Float64},
    amounts::NamedTuple, Gsel::AbstractEffectPrior, Qsel::AbstractEffectPrior,
    X::Matrix{Float64}=ones(amounts.individuals, 1))::NamedTuple

    # Get amounts
    individuals::Int = amounts.individuals
    genes::Int = amounts.genes
    loci::Int = amounts.loci

    # The number of parallel tasks is limited by the number of threads and genes
    nsplit::Int = min(nthreads(), genes)

    # Fixed effects covariates (not centered — used as-is)
    Xcols = [copy(x) for x in eachcol(X)]
    ssqX = [dot(x, x) for x in Xcols]
    fixedeffects = [zeros(size(X, 2)) for _ in 1:genes]

    o = (
        in=(
            Y=mat2vv(Y), # Centered gene expressions
            M=mat2vv(M), # Centered genotypes
            X=Xcols, # Fixed effect covariates (not centered)
            R=[copy(y) for y in eachcol(Y)] # Residuals
        ),
        out=(
            base=(
                fixedeffects=fixedeffects, # Fixed effect coefficients, one vector per gene
                residual_variance=[var(y) for y in eachcol(Y)], # The vector of residual variances
                residual_df=4.0 .* ones(genes), # The vector of degrees of freedom for the prior variance of residuals
                residual_variancescale=[2.0 * var(y) for y in eachcol(Y)] # The vector of prior degrees of freedom times prior residual variances
            ),
            G=_initstate(Gsel, genes, genes; diag=false),
            Q=_initstate(Qsel, genes, loci; diag=true)
        ),
        tmpshared=(
            nsplit=nsplit,
            split_chunks=split_vector(axes(Y, 2), nsplit),
            ssqX=ssqX) # Squared sums of fixed effect covariates
    )

    return o
end

"""

    insert_priors!(d, priors)

Insert priors into the output part of the data bucket. Any field of the
`out` part can be altered: the key for a base field is the field name
(`Rvar`, `Rvardf`, `Rvarscale`), and the key for a prior state field
is the group letter followed by the field name (`Gvar`, `Gvardf`, `Qpi`,
`Glambda2`, ...). Values are broadcast into the fields, so scalars are
recycled across genes. After insertion, the prior variance scales are
recalculated.

"""
function insert_priors!(d::NamedTuple, priors::NamedTuple)::Nothing

    # Insert base priors
    basefields = (Rvar=:residual_variance, Rvardf=:residual_df, Rvarscale=:residual_variancescale)
    for k in keys(priors)
        if haskey(basefields, k)
            getfield(d.out.base, basefields[k]) .= priors[k]
        end
    end

    # Insert prior state priors
    for g in (:G, :Q)
        st = d.out[g]
        for f in fieldnames(typeof(st))
            k = Symbol(g, f)
            if haskey(priors, k)
                getfield(st, f) .= priors[k]
            end
        end
        haskey(priors, Symbol(g, :varscale)) || _varscale!(st)
    end

    # Recalculate the residual variance scale unless it was given as a prior
    haskey(priors, :Rvarscale) || (d.out.base.residual_variancescale .=
        d.out.base.residual_variance .* (d.out.base.residual_df .- 2.0))

    return nothing
end
