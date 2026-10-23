        
        
        
function accuracy(T::AbstractString, E::AbstractString)
        Tmat = readdlm(T, ',', Float64, '\n')
        Emat = reshape(posterior_mean(E), size(Tmat))

        return cor(vec(Tmat),vec(Emat))
end

function posterior_mean(path::AbstractString; delim = '\t')
    
    # Initializing objects
    nlines = 0
    o = Union{Float64,Missing}[]
    
    # Open as gzip if the file ends in .gz, else plain text
    io = endswith(path, ".gz") ? GZip.open(path) : open(path)
    
    try
        # Scrolling the file
        for line in eachline(io)
            isempty(strip(line)) && continue # Skip the line if it is empty
            values = tryparse.(Float64, split(strip(line), delim))
            any(isnothing, values) && error("Non-numeric value on line $(nlines + 1)")
            if nlines == 0
                o = Float64.(values)
            else
                o .+= values
            end
            nlines += 1
        end
    finally
        close(io)
    end
    nlines == 0 && error("No data rows found in $path")
    return o ./ nlines
end