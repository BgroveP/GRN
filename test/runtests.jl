using GRN
using Random
using LinearAlgebra
using Statistics
using Test
using GZip

# Helper: read an output file as a matrix with a given number of values per line
function readoutput(p, ncols)
    isfile(p) || error("Missing output file $p")
    v = [parse(Float64, t) for ln in readlines(p) for t in split(ln, '\t') if !isempty(t)]
    return reshape(v, ncols, :)'
end

# Simulate data with a known network. B[i, j] is the effect of gene i on
# gene j, F[l, j] is the effect of marker l on gene j.
function simulate(n=200, m=5, u=10; seed=1)
    rng = MersenneTwister(seed)
    M = randn(rng, n, u)
    B = zeros(m, m)
    F = zeros(u, m)
    for j in 2:m
        B[j-1, j] = iseven(j) ? 0.8 : -0.8
        j - 1 <= u && (F[j-1, j] = iseven(j) ? 1.0 : -1.0)
    end
    mu = collect(1.0:m)
    Y = zeros(n, m)
    for j in 1:m
        Y[:, j] = mu[j] .+ Y * B[:, j] .+ M * F[:, j] .+ 0.25 .* randn(rng, n)
    end
    return Y, M, B, F
end

@testset "GRN.jl" begin

    @testset "Input checks" begin
        out = joinpath(mktempdir(), "grn")
        @test_throws AssertionError grn(randn(5, 3), randn(6, 3), 2:2:4; outpath=out)
        @test_throws AssertionError grn(randn(10, 2), randn(10, 3), 4:2:2; outpath=out)
        @test_throws AssertionError grn(randn(10, 2), randn(10, 3), 0:2:4; outpath=out)
    end

    @testset "Prior selection" begin
        out = joinpath(mktempdir(), "grn")
        Y = randn(20, 3)
        @test_throws ArgumentError grn(Y, 2:2:4; priors=(Gprior=0.1,), outpath=out)
        @test_throws ArgumentError grn(Y, 2:2:4; priors=(Qprior="spikeslab",), outpath=out)
    end

    @testset "Scalar samplers" begin
        Random.seed!(1)
        # Gamma
        x = [GRN._rand_gamma(4.5) for _ in 1:100000]
        @test isapprox(mean(x), 4.5; rtol=0.02)
        @test isapprox(var(x), 4.5; rtol=0.05)
        # Chisq
        x = [GRN._rand_chisq(7.0) for _ in 1:100000]
        @test isapprox(mean(x), 7.0; rtol=0.02)
        @test isapprox(var(x), 14.0; rtol=0.05)
        # Beta
        x = [GRN._rand_beta(2.0, 5.0) for _ in 1:100000]
        @test isapprox(mean(x), 2.0 / 7.0; rtol=0.02)
        # Inverse Gaussian
        x = [GRN._rand_invgauss(2.0, 3.0) for _ in 1:100000]
        @test isapprox(mean(x), 2.0; rtol=0.05)
    end

    @testset "Priors are inserted" begin
        Y, M, = simulate(50, 3, 4)
        d = GRN.create_databucket(Y, M, (individuals=50, genes=3, loci=4),
            GRN.SpikeSlab(), GRN.SpikeSlab())
        GRN.insert_priors!(d, (mu=[1.0, 2.0, 3.0], Rvar=0.5, Gpi=0.2, Gvar=0.3,
            Gvardf=6.0, Qpi=0.4, Rvardf=5.0))
        @test d.out.base.average == [1.0, 2.0, 3.0]
        @test all(d.out.base.residual_variance .== 0.5)
        @test all(d.out.G.pi .== 0.2)
        @test all(d.out.G.var .== 0.3)
        @test all(d.out.G.varscale .== 0.3 * (6.0 - 2.0))
        @test all(d.out.Q.pi .== 0.4)
        @test all(d.out.base.residual_variancescale .== 0.5 * (5.0 - 2.0))
        # The G prior excludes the diagonal, the Q prior includes all loci
        @test [length(c) for c in d.out.G.candidates] == [2, 2, 2]
        @test [length(c) for c in d.out.Q.candidates] == [4, 4, 4]
    end

    @testset "Analysis and output" begin
        testdir = mktempdir()
        Random.seed!(2)
        Y, M, B, = simulate()

        # Run all prior combinations
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "ss"))
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "lr"),
            priors=(Gprior=Laplace(), Qprior=Ridge()))
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "rl"),
            priors=(Gprior=Ridge(), Qprior=Laplace()))

        # Spike-and-slab output
        @test size(readoutput(joinpath(testdir, "ss", "G_effects.tsv"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "ss", "G_indicator.tsv"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "ss", "G_pi.tsv"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "ss", "G_variance.tsv"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "ss", "Q_effects.tsv"), 50)) == (51, 50)
        @test size(readoutput(joinpath(testdir, "ss", "base_average.tsv"), 5)) == (51, 5)
        @test all(isfinite, readoutput(joinpath(testdir, "ss", "G_variance.tsv"), 5))

        # Laplace output
        @test size(readoutput(joinpath(testdir, "lr", "G_effects.tsv"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "lr", "G_tau2.tsv"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "lr", "G_lambda2.tsv"), 5)) == (51, 5)
        @test all(readoutput(joinpath(testdir, "lr", "G_tau2.tsv"), 25) .> 0)
        @test all(readoutput(joinpath(testdir, "lr", "G_lambda2.tsv"), 5) .> 0)

        # Ridge output
        @test size(readoutput(joinpath(testdir, "rl", "G_effects.tsv"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "rl", "G_variance.tsv"), 5)) == (51, 5)

        # Analysis without genotypes writes no marker-gene files
        grn(Y, 500:10:1000, outpath=joinpath(testdir, "nom"))
        @test isfile(joinpath(testdir, "nom", "G_effects.tsv"))
        @test !isfile(joinpath(testdir, "nom", "Q_effects.tsv"))
        @test size(readoutput(joinpath(testdir, "nom", "G_effects.tsv"), 25)) == (51, 25)

        # The network is recovered with the spike-and-slab prior
        Random.seed!(3)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "rec"),
            priors=(Gpi=0.2, Qpi=0.05))
        Gflat = vec(mean(readoutput(joinpath(testdir, "rec", "G_effects.tsv"), 25), dims=1))
        Ghat = reshape(Gflat, 5, 5)' # Ghat[j, q] is the effect of gene q on gene j
        @test cor(vec(Ghat), vec(B')) > 0.7
    end

    @testset "Allocations" begin
        # The per-gene Gibbs updates must not allocate
        Random.seed!(4)
        Y, M, = simulate(100, 4, 6)
        d = GRN.create_databucket(Y, M, (individuals=100, genes=4, loci=6),
            GRN.SpikeSlab(), GRN.SpikeSlab())
        ssqs = (M=[dot(m, m) for m in d.in.M], Y=[dot(y, y) for y in d.in.Y])
        rv = d.out.base.residual_variance

        f() = begin
            for c in 1:4
                GRN._sample_mean!(d, c, 100)
                GRN._sample_residualvariance!(d, c)
                GRN._sample_effects_gene!(d.out.G, c, d.in.Y, ssqs.Y, d.in.R, rv)
                GRN._sample_hyperparameters_gene!(d.out.G, c, rv)
                GRN._sample_effects_gene!(d.out.Q, c, d.in.M, ssqs.M, d.in.R, rv)
                GRN._sample_hyperparameters_gene!(d.out.Q, c, rv)
            end
        end
        f(); f()
        @test (@allocated f()) == 0

        # One iteration only allocates the tasks spawned by @threads
        f2() = GRN._sample_iteration!(d, ssqs, 100, 6)
        f2(); f2()
        @test (@allocated f2()) < 6 * 4096
    end
end
