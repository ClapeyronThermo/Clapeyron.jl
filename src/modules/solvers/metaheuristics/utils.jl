"""
Shared utilities for metaheuristic solvers.

Only includes helpers currently used by RDEx.
"""

function _sample_index_from_weights(rng::AbstractRNG, weights::AbstractVector{Float64}, len::Int)
    total = 0.0
    @inbounds for i = 1:len
        total += weights[i]
    end
    if !(total > 0.0)
        return 1
    end
    r = rand(rng) * total
    acc = 0.0
    @inbounds for i = 1:len
        acc += weights[i]
        if r <= acc
            return i
        end
    end
    return len
end

@inline function _uniform(rng::AbstractRNG, a::Float64, b::Float64)
    return rand(rng) * (b - a) + a
end
