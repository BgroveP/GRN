
        
        
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
    nlines = 0
    counts = Int[]

    io = endswith(path, ".gz") ? GZip.open(path) : open(path)

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
