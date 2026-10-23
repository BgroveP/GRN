module GRN

# Used packages
using Base.Threads
using LinearAlgebra
using Printf
using Random
using Statistics
using DelimitedFiles
using GZip

# Set number of BLAS threads to one
BLAS.set_num_threads(1)

# Includes
include("postanalysis.jl")
include("priors.jl")
include("initialize.jl")
include("check.jl")
include("estimate.jl")
include("io.jl")
include("sample.jl")
include("utility.jl")

# Exports
export grn
export SpikeSlab, Ridge, Laplace
export accuracy, posterior_mean, posterior_inclusion_probability, derive_indicator

end
