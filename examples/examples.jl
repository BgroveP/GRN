
# Load packages
using DelimitedFiles
using GRN

# Load the input
pkgroot = dirname(dirname(pathof(GRN)))
Y = readdlm("$(pkgroot)/data/Y.csv", ',');
M = readdlm("$(pkgroot)/data/M.csv", ',');

# Spike-and-slab priors, no other priors
grn(Y, M, 10000:50:100000, outpath = "examples/grn_example", overwrite = true)

# Spike-and-slab priors, with hyperparameter priors
grn(Y, M, 10000:50:100000, outpath = "examples/grn_example2", overwrite = true, priors = (Gpi = 0.01, Qpi = 0.01, Rvar = 0.01))

# Bayesian lasso prior for the gene-gene effects, ridge prior for the
# marker-gene effects
grn(Y, M, 10000:50:100000, outpath = "examples/grn_example3", overwrite = true, priors = (Gprior = Laplace(), Qprior = Ridge()))

# Pure Ridge regression
grn(Y, M, 10000:50:100000, outpath = "examples/grn_example4", overwrite = true, priors = (Gprior = Ridge(), Qprior = Ridge()))


# Accuracies
@info "
    Accuracies:
        Gene effects:
            Spike and slab (no priors)       $(round(GRN.accuracy("data/B.csv", "examples/grn_example/G_effects.tsv"), digits = 2))
            Spike and slab (with priors)     $(round(GRN.accuracy("data/B.csv", "examples/grn_example2/G_effects.tsv"), digits = 2))
            Lasso+                           $(round(GRN.accuracy("data/B.csv", "examples/grn_example3/G_effects.tsv"), digits = 2))
            Ridge                            $(round(GRN.accuracy("data/B.csv", "examples/grn_example4/G_effects.tsv"), digits = 2))

        Marker effects:
            Spike and slab (no priors)       $(round(GRN.accuracy("data/F.csv", "examples/grn_example/Q_effects.tsv"), digits = 2))
            Spike and slab (with priors)     $(round(GRN.accuracy("data/F.csv", "examples/grn_example2/Q_effects.tsv"), digits = 2))
            Lasso+                           $(round(GRN.accuracy("data/F.csv", "examples/grn_example3/Q_effects.tsv"), digits = 2))
            Ridge                            $(round(GRN.accuracy("data/F.csv", "examples/grn_example4/Q_effects.tsv"), digits = 2))
            "
