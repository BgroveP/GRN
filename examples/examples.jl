
# Load packages
using DelimitedFiles
using GRN

# Load the input
pkgroot = dirname(dirname(pathof(GRN)))
Y = readdlm("$(pkgroot)/data/Y.csv", ',');
M = readdlm("$(pkgroot)/data/M.csv", ',');

# Spike-and-slab priors, no other priors
grn(Y, M, 10000:50:100000, outpath = "grn_example")

# Spike-and-slab priors, with hyperparameter priors
grn(Y, M, 10000:50:100000, outpath = "grn_example2", priors = (Gpi = 0.01, Qpi = 0.01, Rvar = 0.01))

# Bayesian lasso prior for the gene-gene effects, ridge prior for the
# marker-gene effects
grn(Y, M, 10000:50:100000, outpath = "grn_example3", priors = (Gprior = Laplace(), Qprior = Ridge()))
