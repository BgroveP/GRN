# Functions for creating output files and writing parameters to them

"""

    makeout(p, outputs; force=false)

Create an output folder and open one tab-separated output file per entry in
`outputs`, a NamedTuple that maps a file name to the object that will be
written to it. Returns a vector of NamedTuples with the name, the output
object, and the IO stream.

"""
function makeout(p::AbstractString, outputs::NamedTuple; force::Bool=false, compress::Bool=true)::Vector{NamedTuple}

    # Handle existing directory
    if isdir(p)
        force || error("The output directory $p already exists")
        rm(p, force=true, recursive=true)
        sleep(1) # Wait for the removal to be picked up by the file system
    end

    # Make new directory
    mkpath(p)

    # Create one IO stream per output object
    return [(name=k, values=outputs[k], io= compress ? GZip.open(joinpath(abspath(p), "$k.tsv.gz"), "w", 1) : open(joinpath(abspath(p), "$k.tsv"), "w"))
        for k in keys(outputs)]
end

function _write_row(io::IO, x::AbstractVector{<:Real}, sep::Char)::Nothing
    @inbounds for i in eachindex(x)
        i != firstindex(x) && print(io, sep)
        print(io, x[i])
    end
    print(io, '\n')
    return nothing
end

function _write_row(io::IO, x::AbstractVector{<:AbstractVector}, sep::Char)::Nothing
    first::Bool = true
    for v in x
        for val in v
            first || print(io, sep)
            print(io, val)
            first = false
        end
    end
    print(io, '\n')
    return nothing
end

function _print(i::NamedTuple)::Nothing
    println("Number of individuals:  $(i.individuals)")
    println("Number of genes:        $(i.genes)")
    println("Number of markers:      $(i.loci)")
    println("Number of threads:      $(nthreads())")
    return nothing
end

function close_outputs(outs::Vector{NamedTuple})
    for o in outs
        close(o.io)
    end
end


function write_output(o::NamedTuple; sep::Char='\t')::Nothing
    buf = IOBuffer()
    _write_row!(buf, o.values, sep)
    write(o.io, take!(buf))
    return nothing
end

function _write_row!(buf::IOBuffer, x::AbstractVector{<:Real}, sep::Char)
    for i in eachindex(x)
        i != firstindex(x) && print(buf, sep)
        print(buf, x[i])
    end
    println(buf)
end

function _write_row!(buf::IOBuffer, x::AbstractVector{<:AbstractVector}, sep::Char)
    first = true
    for v in x, val in v
        first || print(buf, sep)
        print(buf, val)
        first = false
    end
    println(buf)
end
