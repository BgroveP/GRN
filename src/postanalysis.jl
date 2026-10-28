
        
# --- Binary file helpers ---

"""
    _is_binary(path)

Return `true` if `path` has a binary extension (`.bin` or `.bin.gz`).
"""
_is_binary(path::AbstractString)::Bool = endswith(path, ".bin") || endswith(path, ".bin.gz")

"""
    _open_output(path)

Open an output file for reading, transparently handling gzip compression.
"""
_open_output(path::AbstractString)::IO = endswith(path, ".gz") ? GZip.open(path) : open(path)

"""
    _read_bin_header(io)

Read and validate the 16-byte binary header from an IO stream.
Returns the number of Float64 values per row (`ncols`).
"""
function _read_bin_header(io::IO)::Int
    magic = read(io, UInt32)
    magic == GRN._BIN_MAGIC || error(
        "Not a GRN binary file (bad magic: 0x$(string(magic, base=16, pad=8)))")
    version = read(io, UInt32)
    version == GRN._BIN_VERSION || error(
        "Unsupported GRN binary version: $version (expected $(GRN._BIN_VERSION))")
    ncols = read(io, Int64)
    ncols > 0 || error("Invalid GRN binary file: ncols = $ncols")
    return ncols
end

"""
    read_binary(path)

Read a GRN binary output file (written with `format=:bin`) and return a
matrix with one row per MCMC iteration and one column per parameter value.
Handles both plain `.bin` and gzip-compressed `.bin.gz` files.

# Example
```julia
# After running grn with format=:bin:
G = GRN.read_binary("output/G_effects.bin")       # 51 × 25 matrix
Q = GRN.read_binary("output/Q_effects.bin.gz")    # 51 × 50 matrix
```
"""
function read_binary(path::AbstractString)
    io = _open_output(path)
    try
        ncols = _read_bin_header(io)
        data_bytes = read(io)
        n_bytes = length(data_bytes)
        n_bytes % sizeof(Float64) == 0 || error(
            "Corrupt binary file: data size ($n_bytes bytes) is not a multiple of 8")
        nvals = n_bytes ÷ sizeof(Float64)
        nvals % ncols == 0 || error(
            "Corrupt binary file: $nvals values is not a multiple of ncols=$ncols")
        nrows = nvals ÷ ncols
        return reshape(reinterpret(Float64, data_bytes), ncols, nrows)'
    finally
        close(io)
    end
end

"""
    _posterior_mean_bin(path)

Streaming posterior mean for binary files. Reads row by row without
loading the entire file into memory.
"""
function _posterior_mean_bin(path::AbstractString)::Vector{Float64}
    io = _open_output(path)
    try
        ncols = _read_bin_header(io)
        sum = zeros(Float64, ncols)
        row = Vector{Float64}(undef, ncols)
        nrows = 0
        while !eof(io)
            read!(io, row)
            @inbounds for i in 1:ncols
                sum[i] += row[i]
            end
            nrows += 1
        end
        nrows == 0 && error("No data rows found in $path")
        return sum ./ nrows
    finally
        close(io)
    end
end

"""
    _pip_bin(path)

Streaming posterior inclusion probability for binary files.
"""
function _pip_bin(path::AbstractString)::Vector{Float64}
    io = _open_output(path)
    try
        ncols = _read_bin_header(io)
        counts = zeros(Int, ncols)
        row = Vector{Float64}(undef, ncols)
        nrows = 0
        while !eof(io)
            read!(io, row)
            @inbounds for i in 1:ncols
                row[i] != 0.0 && (counts[i] += 1)
            end
            nrows += 1
        end
        nrows == 0 && error("No data rows found in $path")
        return counts ./ nrows
    finally
        close(io)
    end
end

# --- Public post-analysis functions (auto-detect binary vs TSV) ---

function accuracy(T::AbstractString, E::AbstractString)
        Tmat = readdlm(T, ',', Float64, '\n')
        Emat = reshape(posterior_mean(E), size(Tmat))

        return cor(vec(Tmat),vec(Emat))
end

function posterior_mean(path::AbstractString; delim = '\t')
    _is_binary(path) && return _posterior_mean_bin(path)
    
    # Initializing objects
    nlines = 0
    o = Union{Float64,Missing}[]
    
    # Open as gzip if the file ends in .gz, else plain text
    io = _open_output(path)
    
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

"""
    posterior_inclusion_probability(path; delim='\\t')

Compute the exact posterior inclusion probability (PIP) for each effect
from an effects file written by `grn` with a spike-and-slab prior. In each
MCMC iteration, an excluded effect is written as exactly `0.0` and an
included effect is non-zero, so the PIP is the proportion of rows where
the value at each position is non-zero.

This gives the exact PIP — identical to averaging the indicator file —
since indicator files are not written by `grn`.

# Example
```julia
# After running grn:
pip = GRN.posterior_inclusion_probability("output/G_effects.tsv.gz")
# Reshape to gene x gene matrix:
G_pip = reshape(pip, n_genes, n_genes)
# Effects included in >50% of iterations:
G_strong = G_pip .> 0.5
```
"""
function posterior_inclusion_probability(path::AbstractString; delim::Char='\t')
    _is_binary(path) && return _pip_bin(path)

    nlines = 0
    counts = Int[]

    io = _open_output(path)

    try
        for line in eachline(io)
            isempty(strip(line)) && continue
            values = tryparse.(Float64, split(strip(line), delim))
            any(isnothing, values) && error("Non-numeric value on line $(nlines + 1)")
            vals = Float64.(values)
            if nlines == 0
                counts = zeros(Int, length(vals))
            end
            @inbounds for i in eachindex(vals)
                vals[i] != 0.0 && (counts[i] += 1)
            end
            nlines += 1
        end
    finally
        close(io)
    end
    nlines == 0 && error("No data rows found in $path")
    return counts ./ nlines
end

"""
    derive_indicator(effects_path; threshold=0.5, delim='\\t')

Derive a binary indicator vector from an effects file written by `grn`
with a spike-and-slab prior. An indicator is `true` wherever the posterior
inclusion probability (PIP) exceeds `threshold` (default: 0.5, i.e. the
effect was included in more than half of the MCMC iterations).

This uses the exact PIP computed from non-zero counts in the effects file
(`posterior_inclusion_probability`), so it gives the same result as
thresholding the averaged indicator file — without the indicator file
being written.

# Example
```julia
# Effects included in >50% of iterations:
ind = GRN.derive_indicator("output/G_effects.tsv.gz")
# Effects included in >80% of iterations:
ind_strong = GRN.derive_indicator("output/G_effects.tsv.gz"; threshold=0.8)
```
"""
function derive_indicator(effects_path::AbstractString; threshold::Float64=0.5, delim::Char='\t')
    pip = posterior_inclusion_probability(effects_path; delim=delim)
    return pip .> threshold
end
