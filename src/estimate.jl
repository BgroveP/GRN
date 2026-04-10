"""

    grn(Y, M, chain; priors=(;), outpath="grn", overwrite=false)

Estimate the gene regulatory network using a matrix with gene expression and a matrix with genotypes.


# Dimensions
 - *n*:           The number of individuals
 - *m*:           The number of genes 
 - *u*:           The number of marker loci


# Mandatory arguments
 - `Y`:             An *n* × *m* matrix with gene expressions                                                           (Float64)
 - `M`:             An *n* × *u* matrix with genotypes                                                                  (Float64)

# Optional arguments
 - `chain`:         The MCMC iterations where parameters are written to the output files                                (UnitRange)
                    For example, giving `10000:50:20000` as input will regard iterations 
                       1:9999 as burn-in, and write parameters for iterations 10000,
                       10050, 10100, ..., 20000    
 - `priors`:        The priors given for the analysis                                                                   (NamedTuple)
 - `outpath`:       The relative path to the folder for output files                                                    (Float64)
 - `overwrite`:     Denotes whether an existing outpath should be overwritten                                           (Bool)
             
 
# Priors
 Priors related to the mean and residual
 - `mu`:            An *m* dimensional vector of average gene expressions                                               (Float64)
 - `Rvar`:          An *m* dimensional vector of residual variances                                                     (Float64)
 - `Rvardf`:        An *m* dimensional vector of degrees of freedoms for the prior residual variance                    (Float64)
 - `Rvarsvale`:     An *m* dimensional vector of ...                                                                    (Float64)

 Priors related to gene-gene effects:
 - `G`:             An *m* × *m* dimensional matrix of effects                                                          (Float64)
 - `Gi`:            An *m* dimensional vector of **m** dimensional vectors of whether an effect is not zero             (Float64)
 - `Gpi`:           An *m* dimensional vector of the probability that an effect is not zero                             (Float64)
 - `Gvar`:          An *m* dimensional vector of variances for effects                                                  (Float64)
 - `Gvardf`:        An *m* dimensional vector of degrees of freedom for the prior variance                              (Float64)
 - `Gvarscale`:     An *m* dimensional vector of ...                                                                    (Float64)

 Priors related to marker-gene effects
 - `Q`:             An *m* × *m* dimensional matrix of effects                                                          (Float64)
 - `Qi`:            An *m* dimensional vector of *m* dimensional vectors of whether an effect is not zero               (Float64)
 - `Qpi`:           An *m* dimensional vector of the probability that an effect is not zero                             (Float64)
 - `Qvar`:          An *m* dimensional vector of variances for effects                                                  (Float64)
 - `Qvardf`:        An *m* dimensional vector of degrees of freedom for the prior variance                              (Float64)
 - `Qvarscale`:     An *m* dimensional vector of ...                                                                    (Float64)

 Priors related to limiting the candidate pathways in the regulatory network
 - `candidates`:    A NamedTuple with two keys (G and Q) with vectors of vectors for indices of candidate genes and markers.
                    Default: (G = repeat([1:*m*], *m*), Q = repeat([1:*u*], *m*))

# Returns
 - `Nothing`: The relevant output is written to files

"""
function grn(Y,
    M,
    chain=100:10:200;
    priors=(;),
    outpath="grn",
    overwrite=false)

    # Get amounts of individuals, genes, and marker loci for later use
    individuals, genes = size(Y)
    loci = size(M, 2)
    println.(
        ["Number of individuals: \t\t$individuals",
         "Number of genes:       \t\t$genes",
         "Number of markers:     \t\t$loci"]
    )

    # Create folder and IO streams
    ## The vector below controls which objects from the databucket that is written to files
    parameternames = ["G", "Gi", "Gpi", "Gvar", "Q", "Qi", "Qpi", "Qvar", "Rvar", "mu"]

    ## This function creates a named collection of IOStreams for later use
    ofiles = makeout(outpath, parameternames, force=overwrite)

    # Checks 
    ## These are not implemented yet, but there should be checks that the number of individuals is the same in M and Y

    # Initialize parameters
    ## This creates a large named tuple with input, output, and temporary arrays (See the io.jl file for details)
    data = create_databucket(Y, M, genes, loci, individuals)

    # Insert priors
    ## This function inserts priors into the output part of the data bucket. 
    ## Any object in the output part of the data bucket can be altered with a prior
    insert_priors!(data, priors)

    # Calculate squared sums for each column of Y and M 
    ssqs = (
        M=dot.(eachcol(data.in.M), eachcol(data.in.M)),
        Y=dot.(eachcol(data.in.Y), eachcol(data.in.Y))
    )

    # Predefine constant distributions
    ##### check julia's randn() function. #### rand!(v) and randn!(v) seem equivalent 
    Normdist = Normal(0, 1) # This is the standard normal distribution. I can just sample from this and scale the output
    Rvaridist = Chisq((data.out.Rvardf...) + individuals) # This is the chisq distribution for residual variances. This is constant across iterations

    # Start analyses
    ## Showprogress gives a progress bar and a timer. I have not experienced slowdowns from it, but could be a culprit.
    @showprogress for i in 1:last(chain)

        # Sample parameters
        _sample_mean!(data, Normdist, individuals)
        _sample_residualvariance!(data, Rvaridist)
        _sample_effectparameters!(data, :G, :Y, ssqs, Normdist)
        _sample_effectparameters!(data, :Q, :M, ssqs, Normdist)

        # Write to file if iteration is an output iteration
        if i in chain
            # Write each output object to its respective file
            @floop for k in keys(ofiles)
                write_output(data, ofiles, k)
            end
        end
    end

    # Close the output files
    for o in ofiles
        close(o)
    end

    # Return
    return nothing
end

