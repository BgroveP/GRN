function create_databucket(Y, M, genes, loci)

    nsplit = min(nthreads(), size(Y, 2))

    return (in=(
            Y=copy(Y) .- mean(Y, dims=1),
            M=copy(M) .- mean(M, dims=1),
            R=copy(Y),
        ),
        out=(
            mu=zeros(genes),
            G=zeros(Float64, genes, genes),
            Gi=zeros(Bool, genes, genes),
            Gpi=0.5*ones(genes),
            Gvar=ones(genes),
            Gvardf=4 * ones(genes),
            Gvarscale=ones(genes),
            Q=zeros(Float64, loci, genes),
            Qi=zeros(Bool, loci, genes),
            Qpi=0.5*ones(genes),
            Qvar=ones(genes),
            Qvardf=4 * ones(genes),
            Qvarscale=ones(genes),
            Rvar=var.(eachcol(Y)),
            Rvardf=4 * ones(genes),
            Rvarscale=ones(genes),
            candidates = (G = repeat([1:genes], genes), Q = repeat([1:loci], genes))
        ),
        tmpshared=(GFloat1=zeros(genes), GFloat2=zeros(genes), GFloat3=zeros(genes), QFloat1=zeros(loci), wlock = ReentrantLock()),
        tmpthread=ntuple(i -> (GFloat1=zeros(genes),
                GFloat2=zeros(genes),
                GFloat3=zeros(genes),
                GFloat4=zeros(genes),
                QFloat1=zeros(loci),
                QFloat2=zeros(loci),
                QFloat3=zeros(loci),
                QFloat4=zeros(loci),
                rng = Random.Xoshiro(rand(1:typemax(Int))),
                values = (dotprod = zeros(1), prob = zeros(1), logD0 = zeros(1), logD1 = zeros(1))), nsplit))
end

function insert_priors!(d::NamedTuple, p::NamedTuple)

    ps = keys(d.out)
    for k in keys(p)
        if in(k, ps)
            data.out[k] .= p[k]
        else
            error("The prior $(k) is not valid")
        end
    end

    return nothing
end


