using GRN
using Random
using LinearAlgebra
using Statistics
using Test
using GZip

# Helper: read an output file as a matrix with a given number of values per line.
# Handles both plain .tsv and gzip-compressed .tsv.gz files.
function readoutput(p, ncols)
    isfile(p) || error("Missing output file $p")
    io = endswith(p, ".gz") ? GZip.open(p) : open(p)
    try
        v = [parse(Float64, t) for ln in eachline(io) for t in split(ln, '\t') if !isempty(t)]
        return reshape(v, ncols, :)'
    finally
        close(io)
    end
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

        # Run all prior combinations (compressed by default)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "ss"))
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "lr"),
            priors=(Gprior=Laplace(), Qprior=Ridge()))
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "rl"),
            priors=(Gprior=Ridge(), Qprior=Laplace()))

        # Spike-and-slab output (compressed .tsv.gz, no indicator file)
        @test size(readoutput(joinpath(testdir, "ss", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test !isfile(joinpath(testdir, "ss", "G_indicator.tsv.gz"))
        @test !isfile(joinpath(testdir, "ss", "G_effects.tsv"))
        @test size(readoutput(joinpath(testdir, "ss", "G_pi.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "ss", "G_variance.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "ss", "Q_effects.tsv.gz"), 50)) == (51, 50)
        @test size(readoutput(joinpath(testdir, "ss", "base_average.tsv.gz"), 5)) == (51, 5)
        @test all(isfinite, readoutput(joinpath(testdir, "ss", "G_variance.tsv.gz"), 5))

        # Laplace output
        @test size(readoutput(joinpath(testdir, "lr", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "lr", "G_tau2.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "lr", "G_lambda2.tsv.gz"), 5)) == (51, 5)
        @test all(readoutput(joinpath(testdir, "lr", "G_tau2.tsv.gz"), 25) .> 0)
        @test all(readoutput(joinpath(testdir, "lr", "G_lambda2.tsv.gz"), 5) .> 0)

        # Ridge output
        @test size(readoutput(joinpath(testdir, "rl", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "rl", "G_variance.tsv.gz"), 5)) == (51, 5)

        # Analysis without genotypes writes no marker-gene files
        grn(Y, 500:10:1000, outpath=joinpath(testdir, "nom"))
        @test isfile(joinpath(testdir, "nom", "G_effects.tsv.gz"))
        @test !isfile(joinpath(testdir, "nom", "Q_effects.tsv.gz"))
        @test size(readoutput(joinpath(testdir, "nom", "G_effects.tsv.gz"), 25)) == (51, 25)

        # The network is recovered with the spike-and-slab prior
        Random.seed!(3)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "rec"),
            priors=(Gpi=0.2, Qpi=0.05))
        Gflat = vec(mean(readoutput(joinpath(testdir, "rec", "G_effects.tsv.gz"), 25), dims=1))
        Ghat = reshape(Gflat, 5, 5)' # Ghat[j, q] is the effect of gene q on gene j
        @test cor(vec(Ghat), vec(B')) > 0.7
    end

    @testset "Compressed output (default)" begin
        testdir = mktempdir()
        Random.seed!(5)
        Y, M, = simulate()

        # Run without specifying compress — default is compressed
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "default"))

        # .tsv.gz files exist, plain .tsv files do not
        @test isfile(joinpath(testdir, "default", "G_effects.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "G_pi.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "G_variance.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "Q_effects.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "base_average.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "base_residual_variance.tsv.gz"))
        @test !isfile(joinpath(testdir, "default", "G_effects.tsv"))
        @test !isfile(joinpath(testdir, "default", "Q_effects.tsv"))

        # Indicator files are never written (PIP derived from effects)
        @test !isfile(joinpath(testdir, "default", "G_indicator.tsv.gz"))
        @test !isfile(joinpath(testdir, "default", "Q_indicator.tsv.gz"))

        # Compressed files can be read back with correct dimensions
        @test size(readoutput(joinpath(testdir, "default", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "default", "G_pi.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "default", "G_variance.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "default", "Q_effects.tsv.gz"), 50)) == (51, 50)
        @test size(readoutput(joinpath(testdir, "default", "base_average.tsv.gz"), 5)) == (51, 5)

        # All values are finite
        @test all(isfinite, readoutput(joinpath(testdir, "default", "G_effects.tsv.gz"), 25))
        @test all(isfinite, readoutput(joinpath(testdir, "default", "G_variance.tsv.gz"), 5))

        # posterior_mean works on .gz files
        pm = GRN.posterior_mean(joinpath(testdir, "default", "G_effects.tsv.gz"))
        @test length(pm) == 25
        @test all(isfinite, pm)

        # Compressed files are smaller than uncompressed
        Random.seed!(5)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "plain"), compress=false)
        gz_size = filesize(joinpath(testdir, "default", "G_effects.tsv.gz"))
        plain_size = filesize(joinpath(testdir, "plain", "G_effects.tsv"))
        @test gz_size < plain_size
    end

    @testset "Plain output (compress=false)" begin
        testdir = mktempdir()
        Random.seed!(6)
        Y, M, = simulate()

        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "plain"), compress=false)

        # Plain .tsv files exist, .tsv.gz files do not
        @test isfile(joinpath(testdir, "plain", "G_effects.tsv"))
        @test isfile(joinpath(testdir, "plain", "G_pi.tsv"))
        @test isfile(joinpath(testdir, "plain", "Q_effects.tsv"))
        @test isfile(joinpath(testdir, "plain", "base_average.tsv"))
        @test !isfile(joinpath(testdir, "plain", "G_effects.tsv.gz"))
        @test !isfile(joinpath(testdir, "plain", "Q_effects.tsv.gz"))

        # Indicator files are never written
        @test !isfile(joinpath(testdir, "plain", "G_indicator.tsv"))

        # Plain files can be read back
        @test size(readoutput(joinpath(testdir, "plain", "G_effects.tsv"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "plain", "G_pi.tsv"), 5)) == (51, 5)
        @test all(isfinite, readoutput(joinpath(testdir, "plain", "G_effects.tsv"), 25))

        # posterior_mean works on plain files
        pm = GRN.posterior_mean(joinpath(testdir, "plain", "G_effects.tsv"))
        @test length(pm) == 25
        @test all(isfinite, pm)
    end

    @testset "PIP and derive_indicator" begin
        testdir = mktempdir()
        Random.seed!(8)
        Y, M, = simulate()

        # Run analysis (compressed by default, indicator files are never written)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "run"))

        effects_path = joinpath(testdir, "run", "G_effects.tsv.gz")

        # --- posterior_inclusion_probability ---

        pip = GRN.posterior_inclusion_probability(effects_path)
        @test pip isa Vector{Float64}
        @test length(pip) == 25  # 5 genes x 5 predictors
        @test all(0.0 .<= pip .<= 1.0)

        # Verify PIP by manual computation: proportion of non-zero rows per column
        effects_mat = readoutput(effects_path, 25)  # 51 rows x 25 cols
        manual_pip = vec(mean(effects_mat .!= 0.0, dims=1))
        @test pip ≈ manual_pip

        # Diagonal always zero (G prior excludes diagonal)
        diag_idx = [1, 7, 13, 19, 25]
        @test all(pip[diag_idx] .== 0.0)

        # --- derive_indicator (PIP-based) ---

        # Default threshold 0.5
        ind = GRN.derive_indicator(effects_path)
        @test ind isa BitVector
        @test length(ind) == 25
        @test all(ind .== (pip .> 0.5))

        # threshold=0.0: included in at least one iteration
        ind_any = GRN.derive_indicator(effects_path; threshold=0.0)
        @test all(ind_any .== (pip .> 0.0))

        # threshold=0.8: high-confidence inclusions
        ind_strong = GRN.derive_indicator(effects_path; threshold=0.8)
        @test all(ind_strong .== (pip .> 0.8))

        # Diagonal always false
        @test all(.!ind[diag_idx])

        # Works on uncompressed files
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "plain"), compress=false)
        pip_plain = GRN.posterior_inclusion_probability(joinpath(testdir, "plain", "G_effects.tsv"))
        @test length(pip_plain) == 25
        ind_plain = GRN.derive_indicator(joinpath(testdir, "plain", "G_effects.tsv"))
        @test ind_plain isa BitVector
        @test all(ind_plain .== (pip_plain .> 0.5))

        # Marker-gene effects
        q_pip = GRN.posterior_inclusion_probability(joinpath(testdir, "run", "Q_effects.tsv.gz"))
        @test length(q_pip) == 50  # 5 genes x 10 markers
        @test all(0.0 .<= q_pip .<= 1.0)
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
