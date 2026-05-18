function get_population!(x::AbstractMatrix, lb::AbstractVector, ub::AbstractVector, rng::AbstractRNG = Random.GLOBAL_RNG)
    dim, npop = size(x)
    length(lb) == dim || throw(DimensionMismatch("lb length must match size(x, 1)"))
    length(ub) == dim || throw(DimensionMismatch("ub length must match size(x, 1)"))
    @inbounds for j = 1:npop, i = 1:dim
        x[i, j] = _uniform(rng, lb[i], ub[i])
    end
    return x
end

get_population(dims::Tuple{Int,Int}, lb::AbstractVector, ub::AbstractVector, rng::AbstractRNG = Random.GLOBAL_RNG) =
    get_population!(Matrix{Float64}(undef, dims...), lb, ub, rng)

function confine!(s::AbstractMatrix{T}, x::AbstractMatrix{T}, lb::AbstractVector{T},
    ub::AbstractVector{T}, factor::T) where {T}
    @inbounds for j in axes(x, 2)
        @simd for i in axes(x, 1)
            tmp = x[i, j] + s[i, j]
            if tmp < lb[i]
                s[i, j] = (lb[i] - x[i, j]) * factor
            elseif tmp > ub[i]
                s[i, j] = (ub[i] - x[i, j]) * factor
            end
        end
    end
    return s
end

@inline function _uniform(rng::AbstractRNG, a::Float64, b::Float64)
    return rand(rng) * (b - a) + a
end
