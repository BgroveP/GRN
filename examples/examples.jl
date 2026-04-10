
# Load packages
using DelimitedFiles
using GRN

# Load the input
pkgroot = dirname(dirname(pathof(GRN)))
Y = readdlm("$(pkgroot)/data/Y.csv", ',');
M = readdlm("$(pkgroot)/data/M.csv", ',');

# No priors
grn(Y, M, 10000:50:100000, outpath = "grn_example")

# With priors
grn(Y, M, 10000:50:100000, outpath = "grn_example2", priors = (Gpi = 0.01, Qpi = 0.01, Rvar = 0.01))

