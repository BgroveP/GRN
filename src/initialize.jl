"""

    create_databucket(Y, M, genes, loci, individuals)

Create the data container with input, parameters, and temporary arrays for inplace operations


# Arguments
 - `Y`:             An *n* × *m* matrix with gene expressions   (Float64)
 - `M`:             An *n* × *u* matrix with genotypes          (Float64)
 - `genes`:         The number of genes                         (Integer)
 - `loci`:          The number of marker loci                   (Integer)
 - `individuals`:   The number of individuals                   (Integer)

# Returns
 - `NamedTuple`: (in = (; ...), out = (; ...), tmpshared = (; ...), tmpthread = (; ...))

"""
function create_databucket(Y, M, genes, loci, individuals)

    # The number of threads used for the analysis is limited by the number of threads available and the number of genes
    nsplit = min(nthreads(), size(Y, 2))

    return (
        # The input is stored as copies because they are altered through centering (Y and M), or in iterations (R)
        in=(
            Y=copy(Y) .- mean(Y, dims=1), # Gene expressions
            M=copy(M) .- mean(M, dims=1), # Genotypes
            R=copy.(eachcol(Y)),          # Residuals
        ),
        # The parameters are stored in out. All objects in out can be changed with the priors input in the grn function.
        out=(
            mu=zeros(genes), # The vector of average effects
            G=copy.(repeat([zeros(genes)], genes)), # The matrix of gene-gene effects
            Gi=copy.(repeat([zeros(Int8, genes)], genes)), # The matrix of whether a gene-gene effect exists
            Gpi=0.5 * ones(genes), # The vector with probabilities of gene-gene effects
            Gvar=ones(genes), # The vector with variances of gene-gene effects
            Gvardf=4 * ones(genes), # The vector of degrees of freedom for the prior variance of gene-gene effects
            Gvarscale=ones(genes), # The ....
            Q=copy.(repeat([zeros(loci)], genes)), # The matrix of marker-gene effects
            Qi=copy.(repeat([zeros(Int8, loci)], genes)), # The matrix of whether marker-gene effects exist
            Qpi=0.5 * ones(genes), # The vector of probabilities for marker-gene effects
            Qvar=ones(genes), # The vector of variances for marker-gene effects
            Qvardf=4 * ones(genes), # The vector of degrees of freedom for the prior variance of marker-gene effects
            Qvarscale=ones(genes), # The ...
            Rvar=var.(eachcol(Y)), # The vector of residual variances 
            Rvardf=4 * ones(genes), # The vector of degrees of freedom for the prior variance of residuals
            Rvarscale=ones(genes), # The ...
            candidates=(G=indicator2indice(copy.(repeat([collect(1:genes)], genes))), Q=repeat([1:loci], genes)) # The candidate sources of effects for G and Q, respectively
        ),
        # The elements below are containers used for inplace operations. 
        # The naming can't be much more specific, because each container is used for multiple types of information.
        # The naming convention is [SIZE][TYPE][NUMBER].
        # For example, GFloat1 is the first number-of-genes dimensional vector of float values
        ## The containers that are shared across threads and not altered under multi-threading
        tmpshared=(
            GFloat1=zeros(genes),
            GFloat2=zeros(genes),
            GFloat3=zeros(genes),
            QFloat1=zeros(loci)),
        ## The containers that are not shared across threads, and are used under multi-threading
        tmpthread=ntuple(i -> (
                GFloat1=zeros(genes),
                GFloat2=zeros(genes),
                GFloat3=zeros(genes),
                GFloat4=zeros(genes),
                QFloat1=zeros(loci),
                QFloat2=zeros(loci),
                QFloat3=zeros(loci),
                QFloat4=zeros(loci),
                IFloat1=zeros(individuals),
                values=(
                    dotprod=zeros(1),
                    prob=zeros(1),
                    logD0=zeros(1),
                    logD1=zeros(1))), # These are values that are reused across iterations 
            nsplit))
end

function insert_priors!(d::NamedTuple, p::NamedTuple)

    # Insert priors
    ps = keys(d.out)
    for k in keys(p)
        if in(k, ps)
            d.out[k] .= p[k]
        else
            error("The prior $(k) is not valid")
        end
    end

    # Update the prior times degrees of freedom for variances
    d.out.Gvarscale .= d.out.Gvar .* (d.out.Gvardf .- 2.0)
    d.out.Qvarscale .= d.out.Qvar .* (d.out.Qvardf .- 2.0)
    d.out.Rvarscale .= d.out.Rvar .* (d.out.Rvardf .- 2.0)

    return nothing
end


