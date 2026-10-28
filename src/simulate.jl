# Simulation of gene expression data from a gene regulatory network.

"""
    simulate(; n=200, m=10, u=0, seed=nothing, density=0.3, effect_sd=0.5,
             noise_sd=0.5, G=nothing, Q=nothing, mu=nothing, M=nothing,
             X=nothing, Beta=nothing)

Simulate gene expression data and genotypes from a linear gene regulatory
network model.

# Model

Without fixed effects:

    Y[:, j] = mu[j] + Y * G[:, j] + M * Q[:, j] + epsilon

With fixed effects (`X` given):

    Y[:, j] = X * Beta[:, j] + Y * G[:, j] + M * Q[:, j] + epsilon

where `G` is an *m* × *m* strictly upper-triangular matrix of gene-gene
effects (a directed acyclic graph), `Q` is a *u* × *m* matrix of
marker-gene effects, `Beta` is a *p* × *m* matrix of fixed effect
coefficients, and `epsilon` is Gaussian noise with standard deviation
`noise_sd`. `G[i, j]` is the effect of gene `i` on gene `j`.

# Keyword arguments

- `n`:          Number of individuals (default: 200)
- `m`:          Number of genes (default: 10)
- `u`:          Number of marker loci (default: 0, no genotypes)
- `seed`:       Random seed (default: random)
- `density`:    Probability of a non-zero effect (default: 0.3)
- `effect_sd`:  Standard deviation of non-zero effect sizes (default: 0.5)
- `noise_sd`:   Standard deviation of the residual noise (default: 0.5)
- `G`:          Optional *m* × *m* matrix of true gene-gene effects.
                Must be strictly upper triangular (G[i,j] = 0 for i >= j).
                If given, `density` and `effect_sd` are ignored for G.
- `Q`:          Optional *u* × *m* matrix of true marker-gene effects.
                If given, `density` and `effect_sd` are ignored for Q.
- `mu`:         Optional *m*-vector of intercepts (default: `collect(1.0:m)`).
                Ignored when `X` is given.
- `M`:          Optional *n* × *u* matrix of genotypes. If given, `u` is
                inferred from its size.
- `X`:          Optional *n* × *p* matrix of fixed effect covariates.
                When given, the intercept is replaced by `X * Beta[:, j]`.
                Include a column of ones for an intercept.
- `Beta`:       Optional *p* × *m* matrix of true fixed effect coefficients.
                If `X` is given but `Beta` is not, a sparse random matrix is
                generated using `density` and `effect_sd`.

# Returns

- `NamedTuple`: `(Y, M, G, Q, mu, X, Beta)` — the simulated data and true
  parameters. `X` is `zeros(n, 0)` and `Beta` is `zeros(0, m)` when no
  fixed effects are used.

# Example

```julia
# Simulate data from a random network with 200 individuals, 20 genes,
# and 50 marker loci
sim = simulate(n=200, m=20, u=50, seed=42)
grn(sim.Y, sim.M, 1000:100:5000, outpath="results")

# Compare inferred effects with truth
Ghat = reshape(GRN.posterior_mean("results/G_effects.tsv.gz"), 20, 20)
cor(vec(Ghat), vec(sim.G'))

# Simulate with fixed effects (e.g. sex and age)
X = [ones(200) rand(0:1, 200) randn(200)]
sim = simulate(n=200, m=20, u=50, X=X, seed=42)
grn(sim.Y, sim.M, 1000:100:5000; X=X, outpath="results_fe")
```
"""
function simulate(; n::Int=200, m::Int=10, u::Int=0,
    seed::Union{Int, Nothing}=nothing,
    density::Float64=0.3, effect_sd::Float64=0.5, noise_sd::Float64=0.5,
    G::Union{Matrix{Float64}, Nothing}=nothing,
    Q::Union{Matrix{Float64}, Nothing}=nothing,
    mu::Union{Vector{Float64}, Nothing}=nothing,
    M::Union{Matrix{Float64}, Nothing}=nothing,
    X::Union{Matrix{Float64}, Nothing}=nothing,
    Beta::Union{Matrix{Float64}, Nothing}=nothing)::NamedTuple

    rng = seed === nothing ? MersenneTwister() : MersenneTwister(seed)

    # Genotypes: use provided M or generate
    if M !== nothing
        size(M, 1) == n || throw(ArgumentError(
            "M must have $n rows (got $(size(M, 1)))"))
        u = size(M, 2)
    else
        M = u > 0 ? randn(rng, n, u) : zeros(n, 0)
    end

    # Fixed effects covariates
    if X !== nothing
        size(X, 1) == n || throw(ArgumentError(
            "X must have $n rows (got $(size(X, 1)))"))
    end
    p = X !== nothing ? size(X, 2) : 0

    # True gene-gene effects (strictly upper triangular DAG)
    if G === nothing
        G = zeros(m, m)
        for j in 2:m, i in 1:j-1
            if rand(rng) < density
                G[i, j] = effect_sd * randn(rng)
            end
        end
    else
        size(G) == (m, m) || throw(ArgumentError(
            "G must be $m × $m (got $(size(G)))"))
        for j in 1:m, i in j:m
            G[i, j] == 0.0 || throw(ArgumentError(
                "G must be strictly upper triangular (G[$i, $j] must be 0 for i >= j)"))
        end
    end

    # True marker-gene effects
    if Q === nothing
        Q = u > 0 ? (rand(rng, u, m) .< density) .* (effect_sd .* randn(rng, u, m)) :
            zeros(0, m)
    else
        size(Q) == (u, m) || throw(ArgumentError(
            "Q must be $u × $m (got $(size(Q)))"))
    end

    # True fixed effect coefficients
    if X !== nothing
        if Beta === nothing
            Beta = p > 0 ? (rand(rng, p, m) .< density) .* (effect_sd .* randn(rng, p, m)) :
                zeros(0, m)
        else
            size(Beta) == (p, m) || throw(ArgumentError(
                "Beta must be $p × $m (got $(size(Beta)))"))
        end
        mu = zeros(m) # intercept is replaced by X * Beta
    else
        Beta = zeros(0, m)
        if mu === nothing
            mu = collect(1.0:m)
        else
            length(mu) == m || throw(ArgumentError(
                "mu must have length $m (got $(length(mu)))"))
        end
    end

    # Simulate gene expressions column by column (topological order)
    Y = zeros(n, m)
    for j in 1:m
        yj = noise_sd .* randn(rng, n)
        if X !== nothing
            yj .+= X * view(Beta, :, j)
        else
            yj .+= mu[j]
        end
        yj .+= Y * G[:, j]
        u > 0 && (yj .+= M * view(Q, :, j))
        Y[:, j] = yj
    end

    X_out = X !== nothing ? X : zeros(n, 0)
    return (Y=Y, M=M, G=G, Q=Q, mu=mu, X=X_out, Beta=Beta)
end
