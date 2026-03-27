
using DelimitedFiles 
using GRN
using BenchmarkTools
using ProfileView
#using ProgressMeter
include("exclude/GRN_BWGR.jl")

Y = readdlm("data/Y.csv", ',')
M = readdlm("data/M.csv", ',')

maxiter = 10000
burnin = 1000
interleave = 10
grn(Y, M, burnin:interleave:maxiter, overwrite = true, outpath = "grn")

@time estGRN_Gibbs(Y,M,maxiter,burnin,interleave,1.0,1.0;outFolder="grn_control")
@profview grn(Y, M, 1000:10:2000, overwrite = true, outpath = "grn")


