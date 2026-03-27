
function makeout(p::String, fs::Vector{String}; force::Bool=false)

    # Handle existing directory
    if isdir(p) & force
        rm(p, force=force, recursive=true)
        sleep(1)
    elseif isdir(p) & ~force
        error("The output directory folder $p already exists.")
    end

    # Make new directory
    mkpath(p)
    println("Creating the output dir for MCMC samples:   $(abspath(p))")

    # Create output tuple

    return NamedTuple{Tuple(Symbol.(fs))}(Tuple(open.(joinpath.(abspath(p), fs), "w")))
end


function write_output(d, o, i; sep='\t')
    x = d.out[i]
    D = length(axes(x))
    b = o[i]

    for col in axes(x, 2)
        for row in axes(x, 1)
            if D == 2
                print(b, x[row, col], sep)
            else
                print(b, x[row], sep)
            end
    end
end
    print(b, '\n')
    return nothing
end



