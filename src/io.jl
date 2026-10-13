# Functions for creating output files and writing parameters to them

"""

    makeout(p, outputs; force=false)

Create an output folder and open one tab-separated output file per entry in
`outputs`, a NamedTuple that maps a file name to the object that will be
written to it. Returns a vector of NamedTuples with the name, the output
object, and the IO stream.

"""
function makeout(p::AbstractString, outputs::NamedTuple; force::Bool=false)::Vector{NamedTuple}

    # Handle existing directory
    if isdir(p)
        force || error("The output directory $p already exists")
        rm(p, force=true, recursive=true)
        sleep(1) # Wait for the removal to be picked up by the file system
    end

    # Make new directory
    mkpath(p)

    # Create one IO stream per output object
    return [(name=k, values=outputs[k], io=open(joinpath(abspath(p), "$k.tsv"), "w"))
        for k in keys(outputs)]
end

"""

    write_output(o; sep='\\t')

Write one row of the output object in `o` to its IO stream. Vectors are
written as one row of tab-separated values, and vectors of vectors are
written as one flattened row per gene.

"""
function write_output(o::NamedTuple; sep::Char='\t')::Nothing
    _write_row(o.io, o.values, sep)
    return nothing
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
