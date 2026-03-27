function grn(Y,
    M,
    chain=100:10:200;
    priors=(;),
    outpath="grn4",
    overwrite=false)

    # Get amounts 
    individuals, genes = size(Y)
    loci = size(M, 2)
    println.(
        ["Number of individuals: \t\t$individuals",
        "Number of genes:       \t\t$genes",
        "Number of markers:     \t\t$loci"]
    )

    # Create folder and IO streams
    parameternames = ["G", "Gi", "Gpi", "Gvar", "Q", "Qi", "Qpi", "Qvar", "Rvar", "mu"]
    ofiles = makeout(outpath, parameternames, force=overwrite)

    # Checks

    # Initialize parameters
    data = create_databucket(Y, M, genes, loci)

    ## Insert priors
    insert_priors!(data, priors)

    # Calculate squared sums for each column of Y and M 
    ssqs = (
        M = dot.(eachcol(data.in.M), eachcol(data.in.M)),
        Y = dot.(eachcol(data.in.Y), eachcol(data.in.Y))
    )

    # Predefine constant distributions
    Normdist = Normal(0, 1)
    Rvaridist = Chisq((data.out.Rvardf...) + individuals)

    # Calculate scales from priors 
    # THIS IS NOT SCALE BUT df*Scale, Seems correctly used at other places though.
    data.out.Gvarscale .= data.out.Gvar .* (data.out.Gvardf .- 2.0)
    data.out.Qvarscale .= data.out.Qvar .* (data.out.Qvardf .- 2.0)
    data.out.Rvarscale .= data.out.Rvar .* (data.out.Rvardf .- 2.0)

    #Start analyses
    @showprogress for i in 1:last(chain)

        _sample_mean!(data, Normdist, individuals)
        _sample_residualvariance!(data, Rvaridist)
        _sample_effects!(data, :G, :Y, ssqs, Normdist, omitdiag=true)
        _sample_effects!(data, :Q, :M, ssqs, Normdist)

        # Maybe write to file
        if i in chain
            @threads for k in keys(ofiles)
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

