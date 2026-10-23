# Functions for creating output files and writing parameters to them

"""
    makeout(p, outputs; force=false, compress=true)

Create an output folder and open one tab-separated output file per entry in
`outputs`, a NamedTuple that maps a file name to the object that will be
written to it. Returns a vector of NamedTuples with the name, the output
object, the IO stream, and a persistent IOBuffer. Files are gzip-compressed
(level 1) by default, producing `.tsv.gz` files; set `compress=false` for
plain `.tsv` files.
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

    # Create one IO stream and one persistent buffer per output object
    return [(name=k, values=outputs[k],
             io=compress ? GZip.open(joinpath(abspath(p), "$k.tsv.gz"), "w", 1)
                          : open(joinpath(abspath(p), "$k.tsv"), "w"),
             buf=IOBuffer(sizehint=_BUFSIZE))
            for k in keys(outputs)]
end

"""
    write_output(o; sep='\t')

Write one row of the output object in `o` to its IO stream. Vectors are
written as one row of tab-separated values, and vectors of vectors are
written as one flattened row. Uses a persistent IOBuffer (1 MB, created in
`makeout`) that is flushed to the stream whenever it fills, so memory usage
stays bounded even for very large rows (e.g. 17K × 17K effects). `@printf`
is used for fast Float64 formatting.
"""
function write_output(o::NamedTuple; sep::Char='\t')::Nothing
    buf = o.buf
    _write_row!(buf, o.io, o.values, sep)
    write(o.io, take!(buf))  # flush any remaining data
    return nothing
end

# Flush the buffer to the stream if it has reached the size threshold.
@inline function _flush_if_full!(buf::IOBuffer, io::IO)::Nothing
    if position(buf) >= _BUFSIZE
        write(io, take!(buf))
    end
    return nothing
end

# --- Fast path: 1D Float64 vector ---

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{Float64}, sep::Char)::Nothing
    @inbounds for i in eachindex(x)
        i != firstindex(x) && print(buf, sep)
        @printf(buf, "%.10g", x[i])
        _flush_if_full!(buf, io)
    end
    println(buf)
    return nothing
end

# --- Fast path: nested Float64 vectors (flattened into one row) ---

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{<:AbstractVector{Float64}}, sep::Char)::Nothing
    first = true
    for v in x
        @inbounds for val in v
            if first
                @printf(buf, "%.10g", val)
                first = false
            else
                print(buf, sep)
                @printf(buf, "%.10g", val)
            end
        end
        _flush_if_full!(buf, io)
    end
    println(buf)
    return nothing
end

# --- Generic fallbacks for non-Float64 element types ---

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{<:Real}, sep::Char)::Nothing
    first = true
    for val in x
        first || print(buf, sep)
        print(buf, val)
        first = false
        _flush_if_full!(buf, io)
    end
    println(buf)
    return nothing
end

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{<:AbstractVector}, sep::Char)::Nothing
    first = true
    for v in x
        for val in v
            first || print(buf, sep)
            print(buf, val)
            first = false
            _flush_if_full!(buf, io)
        end
    end
    println(buf)
    return nothing
end

function _print(i::NamedTuple)::Nothing
    println("Number of individuals:  $(i.individuals)")
    println("Number of genes:        $(i.genes)")
    println("Number of markers:      $(i.loci)")
    println("Number of threads:      $(nthreads())")
    return nothing
end

"""
    close_outputs(outs)

Close all IO streams in the vector of output NamedTuples. Must be called
to flush gzip trailers and release file handles.
"""
function close_outputs(outs::Vector{NamedTuple})::Nothing
    for o in outs
        close(o.io)
    end
    return nothing
end
