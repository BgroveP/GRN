# Checks for the input to the grn function
# 1: Do Y, M, and X fit in dimensions?
# 2: Is the chain an increasing range?

function _checkinput(Y::Matrix{Float64}, M::Matrix{Float64}, chain::AbstractRange{Int},
    X::Matrix{Float64})::Nothing
    @assert size(Y, 1) >= 2 "Y must have at least two individuals"
    @assert size(Y, 2) >= 2 "Y must have at least two genes"
    if size(M, 2) > 0
        @assert size(M, 1) == size(Y, 1) "The number of rows in Y and M must match"
        @assert all(isfinite, M) "M must only contain finite values"
    end
    @assert size(X, 1) == size(Y, 1) "The number of rows in X and Y must match"
    @assert all(isfinite, X) "X must only contain finite values"
    for (k, x) in enumerate(eachcol(X))
        @assert dot(x, x) > 0.0 "Column $k of X is all zeros"
    end
    @assert step(chain) > 0 "The chain must be an increasing range"
    @assert !isempty(chain) "The chain must contain at least one iteration"
    @assert first(chain) >= 1 "The chain must start at iteration 1 or later"
    @assert all(isfinite, Y) "Y must only contain finite values"
    return nothing
end
