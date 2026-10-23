"""

    grn(Y, M, chain; priors=(;), outpath="grn", overwrite=false)

Estimate the gene regulatory network using a matrix with gene expression and
a matrix with genotypes.


# Dimensions
 - *n*:           The number of individuals
 - *m*:           The number of genes
 - *u*:           The number of marker loci


# Mandatory arguments
 - `Y`:             An *n* × *m* matrix with gene expressions                                                           (Float64)
 - `M`:             An *n* × *u* matrix with genotypes                                                                  (Float64)

# Optional arguments
 - `chain`:         The MCMC iterations where parameters are written to the output files                                (AbstractRange{Int})
                    For example, giving `10000:50:20000` as input will regard iterations
                       1:9999 as burn-in, and write parameters for iterations 10000,
                       10050, 10100, ..., 20000
 - `priors`:        The priors given for the analysis                                                                   (NamedTuple)
 - `outpath`:       The relative path to the folder for output files                                                    (String)
 - `overwrite`:     Denotes whether an existing outpath should be overwritten                                           (Bool)
 - `compress`:      Write gzip-compressed `.tsv.gz` files (default, level 1, fast). Set to `false` for plain `.tsv`.  (Bool)

# Priors
 The prior for gene-gene (G) and marker-gene (Q) effects is selected with
 the `Gprior` and `Qprior` keys, which take a prior selector:

  - `SpikeSlab(; includeprob, var, vardf)`: Spike-and-slab prior with sampled
    inclusion probabilities (default)
  - `Ridge(; var, vardf)`: Gaussian prior on all effects with a sampled
    per-gene prior variance
  - `Laplace(; lambda2, a, b)`: Bayesian lasso prior with per-effect local
    variances and a sampled per-gene penalty `lambda2`

 For example, `priors = (Gprior = Laplace(), Qprior = Ridge())`. All other
 keys set hyperparameters of the model, and are inserted into the state
 before the analysis. The key for a base field is the field name, and the
 key for a prior state field is the group letter followed by the field name.

 Priors related to the mean and residual
 - `mu`:            An *m* dimensional vector of average gene expressions                                               (Float64)
 - `Rvar`:          An *m* dimensional vector of residual variances                                                     (Float64)
 - `Rvardf`:        An *m* dimensional vector of degrees of freedoms for the prior residual variance                    (Float64)
 - `Rvarscale`:     An *m* dimensional vector of prior degrees of freedom times prior residual variances                 (Float64)

 Priors related to gene-gene effects (replace `G` with `Q` for marker-gene effects)
 - `Geffects`:      An *m* × *m* dimensional matrix of effects                                                          (Float64)
 - `Gpi`:          An *m* dimensional vector of the probability that an effect is not zero                             (Float64)
 - `Gvar`:          An *m* dimensional vector of variances for effects                                                  (Float64)
 - `Gvardf`:        An *m* dimensional vector of degrees of freedom for the prior variance                              (Float64)
 - `Gvarscale`:     An *m* dimensional vector of prior degrees of freedom times prior effect variances                  (Float64)
 - `Gtau2`:         An *m* dimensional vector of **m** dimensional vectors of local variances (Laplace)                 (Float64)
 - `Glambda2`:      An *m* dimensional vector of squared penalties (Laplace)                                            (Float64)
 - `Ga`:            An *m* dimensional vector of shapes for the prior of lambda2 (Laplace)                               (Float64)
 - `Gb`:            An *m* dimensional vector of rates for the prior of lambda2 (Laplace)                                (Float64)
 - `Gcandidates`:   An *m* dimensional vector of vectors with indices of candidate genes                                (Int)

# Returns
 - `Nothing`: The relevant output is written to files

"""
function grn(Y::Matrix{Float64},
    M::Matrix{Float64}=zeros(0, 0),
    chain::AbstractRange{Int}=100:10:200;
    priors::NamedTuple=(;),
    outpath::AbstractString="grn",
    overwrite::Bool=false,
    compress::Bool=true)

    # Check input
    _checkinput(Y, M, chain)

    # Get amounts of individuals, genes, and marker loci for later use
    amounts = (individuals=size(Y, 1), genes=size(Y, 2), loci=size(M, 2))

    # Print
    _print(amounts)

    # Select the priors for gene-gene and marker-gene effects
    Gsel = get(priors, :Gprior, SpikeSlab())
    Qsel = get(priors, :Qprior, SpikeSlab())
    Gsel isa AbstractEffectPrior || throw(ArgumentError(
        "Gprior must be a prior selector (SpikeSlab(), Ridge(), or Laplace()), got $(typeof(Gsel))"))
    Qsel isa AbstractEffectPrior || throw(ArgumentError(
        "Qprior must be a prior selector (SpikeSlab(), Ridge(), or Laplace()), got $(typeof(Qsel))"))

    # Initialize parameters
    ## This creates a large named tuple with input, output, and temporary arrays
    data = create_databucket(Y, M, amounts, Gsel, Qsel)

    # Insert priors
    ## This function inserts priors into the output part of the data bucket
    insert_priors!(data, priors)

    # Calculate squared sums for each column of Y and M
    ssqs = (
        M=[dot(m, m) for m in data.in.M],
        Y=[dot(y, y) for y in data.in.Y]
    )

    # Create folder and IO streams
    ## This function creates a named collection of IOStreams for later use
    ofiles = makeout(outpath, _outputs(data, amounts.loci); force=overwrite, compress=compress)

    try
        # Start analysis
        for i in 1:last(chain)

            # Sample parameters
            ## The full Gibbs sweep is performed per gene, in parallel
            ## over contiguous chunks of genes
            _sample_iteration!(data, ssqs, amounts.individuals, amounts.loci)

            # Write to file if iteration is an output iteration
            if i in chain
                @threads for k in eachindex(ofiles)
                    write_output(ofiles[k])
                end
            end
        end
    finally
        # Close the output files
        for o in ofiles
            close(o.io)
        end
    end

    # Return
    return nothing
end

# Convenience method for analyses without genotypes
function grn(Y::Matrix{Float64}, chain::AbstractRange{Int};
    priors::NamedTuple=(;),
    outpath::AbstractString="grn",
    overwrite::Bool=false,
    compress::Bool=true)
    return grn(Y, zeros(0, 0), chain; priors=priors, outpath=outpath, overwrite=overwrite,
        compress=compress)
end

"""

    _outputs(d, loci)

Collect the output objects of the data bucket into a flat NamedTuple with
one entry per output file. The keys become the file names. The marker-gene
group is skipped when there are no marker loci. Indicator matrices are not
included — the posterior inclusion probability can be derived exactly from
the effects files with `GRN.posterior_inclusion_probability`.

"""
function _outputs(d::NamedTuple, loci::Int)::NamedTuple

    # Base outputs
    o = (
        base_average=d.out.base.average,
        base_residual_variance=d.out.base.residual_variance
    )

    # Prior state outputs
    groups = loci > 0 ? (:G, :Q) : (:G,)
    for g in groups
        st = d.out[g]
        so = _outputs(st)
        for k in keys(so)
            o = merge(o, NamedTuple{(Symbol(g, "_", k),)}((so[k],)))
        end
    end

    return o
end
