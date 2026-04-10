function split_vector(x::Union{UnitRange{Int},AbstractVector,Base.OneTo{Int64}}, n::Int)

    otype = typeof(x) == Base.OneTo{Int64} ? UnitRange{Int} : typeof(x)
    y = Vector{otype}(undef, n)
    z::Int = length(x)
    r::Int = z
    for i in eachindex(y)
        thistake = ((z-r+1):z)[1:(Int(ceil(r / (n + 1 - i))))]
        y[i] = x[thistake]
        r = r - length(thistake)
    end
    return y
end


function dot!(z, x, y)
    z .= x .* y
    return sum(z)
end

function indicator2indice(v)
    for i in axes(v, 1)
        b = v[i]
        deleteat!(b, i)
    end
    return v
end

