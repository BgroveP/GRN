        
        
        
function accuracy(T::AbstractString, E::AbstractString)
        Tmat = readdlm(T, ',', Float64, '\n')
        Emat = reshape(mean.(eachcol(readdlm(E, '\t', Float64, '\n'))), size(Tmat))

        return cor(vec(Tmat),vec(Emat))
end
