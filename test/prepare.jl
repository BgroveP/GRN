
# 
# Import packages
using RCall
using DelimitedFiles
using DataFrames
using CSV
using Distributions
using Random
using LinearAlgebra
using Random

# Initialize values
Random.seed!(parse(Int, "123"));
totalGene = 10
candidatePerGene = 20
eQTLPerGene = 4 #perGene
eQTLProp = eQTLPerGene/candidatePerGene #prop in all eQTLs
eQTLVar = 0.01 #for sseMQ to work. But i dont use its simulated values.
repl=1
#intercept is zero, as i set.

# 
include("exclude/GRN_BWGR.jl")
path2Omics = "exclude/mySimData/"
path2Analysis = "exclude/myAnalysis/"

workFolderHOL    = path2Analysis*"HOL/repl$repl"
workFolderJER    = path2Analysis*"JER/repl$repl"
workFolderHOLJER = path2Analysis*"HOLJER/repl$repl"
workFolderJERHOL = path2Analysis*"JERHOL/repl$repl"

eQTLHOL = CSV.read(path2Omics*"/HOL/repl$repl/eQTLHOL", Tables.matrix, header = false)
HOL_IDs = collect(1:size(eQTLHOL, 1))
eQTLHOL = eQTLHOL'

eQTLJER = CSV.read(path2Omics*"/JER/repl$repl/eQTLJER", Tables.matrix, header = false)
JER_IDs = collect(1:size(eQTLJER, 1))
eQTLJER = eQTLJER'

function simulateData(inputData, nGene, nTrueEQTL, propEQTL, varEQTL)
	simOut = rcopy(R"""
			library(ssemQr)
			library(network)
			library(ggnetwork)
			library(igraph)
			library(Matrix)
			library(fssemR)
		N = 1000
		Ng = $nGene                                                          # gene number
		Nk = $nGene * $eQTLPerGene                                           # eQTL number
		Ns = 12/Ng   #12/Ng                                                         # sparsity of GRN
		sigma2 = $varEQTL                                                        # sigma2
		Es = $propEQTL                                                             # sparsity of {\it cis}-eQTL
		set.seed(123)  ####Fixes network!
		data = randomeQTLdata(n = N, p = Ng, k = Nk, sparse =  Ns, sqtl = Es, intercept = 0, sigma2 = sigma2, esize = c(0.5, 5), coefs = c(0.2, 0.6), type = "DG", dag = TRUE, overlap = "none", span = TRUE, noqtl = FALSE, rmv = 0.0)
		print(str(data))
		rownames(data$Vars$B) = colnames(data$Vars$B) = rownames(data$Vars$F) = rownames(data$Data$Y)
		colnames(data$Vars$F) = rownames(data$Data$X)
		GE = get.edgelist(graph.adjacency(t(data$Vars$B) != 0))
		QE = which(t(data$Vars$F) != 0, arr.ind = TRUE)
		QE[,2] = rownames(data$Vars$F)[QE[,2]]
		QE[,1] = rownames(QE)
		GRN = network(rbind(GE, QE), matrix.type = "edgelist", directed = TRUE)
		myB = as.matrix(data$Vars$B,dim(data$Vars$B),nrow=T)
		myF = as.matrix(data$Vars$F,dim(data$Vars$F),nrow=T)
		myMu = as.matrix(data$Vars$Mu,dim(data$Vars$Mu),nrow=T) #not used anywhere
		myX = $inputData #I use my own data
		myE = matrix(rnorm(dim(myX)[2]*Ng,0,sqrt(sigma2)),Ng,dim(myX)[2])
		myY = solve(diag(1,dim(myB))-myB)%*%((myF%*%myX)+myE)
		mySk = data$Data$Sk
		myList = list(myB,myF,myY,myMu,myX,myE,mySk,rownames(data$Data$Y),rownames(data$Data$X))
		names(myList) = c("matrixB","matrixF","myY","mu","myX","myE","mySk","Ynames","Xnames")
		myList
		""")
	return simOut
end

simHOL = simulateData(eQTLHOL, totalGene, eQTLPerGene, eQTLProp, eQTLVar)
simJER = simulateData(eQTLJER, totalGene, eQTLPerGene, eQTLProp, eQTLVar)
candidatePerGene = [(i*20-20+1):(i*20) for i in 1:10]
Ftrue = Matrix{Float64}(simJER[:matrixF]')
Btrue = Matrix{Float64}(simJER[:matrixB]')

M = simJER[:myX]' |> Matrix{Float64}
Y = simJER[:myY]'
chain = 100:10:200
priors = (;)
priors = (Gvar = 0.5, Qvar = 0.5, Rvar = 0.01, Gpi = 0.2, Qpi = 0.2)
candidates = [(i*20-20+1):(i*20) for i in 1:10]
Random.seed!(123)


writedlm("data/Y.csv", Y, ",")
writedlm("data/M.csv", M, ",")