function split_vector(x::AbstractVector{Int}, n::Int)::Vector{Vector{Int}}

    y = Vector{Vector{Int}}(undef, n)
    z::Int = length(x)
    r::Int = z
    for i in eachindex(y)
        thistake = ((z - r + 1):z)[1:(Int(ceil(r / (n + 1 - i))))]
        y[i] = collect(x[thistake])
        r = r - length(thistake)
    end
    return y
end

function mat2vv(X::Matrix{Float64}, center::Bool=true)::Vector{Vector{Float64}}

    vv = [copy(v) for v in eachcol(X)]

    if center
        for i in axes(vv, 1)
            vv[i] .-= mean(vv[i])
        end
    end
    return vv
end
