## This is a vendored version of the code in EoSSuperancillaries.jl v1.4

struct ChebyshevRange{R,T}
    range::R
    coeffs::T
end

const ChebyshevRangeV64 = ChebyshevRange{Vector{Float64},Vector{Vector{Float64}}}
const ChebyshevRangeVec{T} = ChebyshevRange{Vector{T},Vector{Vector{T}}} where T

#overload of searchsortedfirst to also match the first element.
_searchsortedfirst(xx::Tuple,x) = _searchsortedfirst(SVector(xx),x)
_searchsortedfirst(xx,x) = Base.searchsortedfirst(xx,x) + isequal(x,first(xx))

#evaluation of ranges of chebyshev coefficients
function cheb_eval(Base.@specialize(cheb), Base.@specialize(x̃))
    Cₙ = cheb.coeffs
    x̃range = cheb.range
    x̃min = first(x̃range)
    x̃max = last(x̃range)
    if !(x̃min <= x̃ <= x̃max)
        R = promote_type(eltype(x̃range), eltype(eltype(Cₙ)),typeof(x̃))
        #x is not in range
        return zero(R)/zero(R)
    end
    x̄,i = cheb_xrange(x̃range,x̃)
    Cₙi = Cₙ[i]
    return cheb_eval(Cₙi,x̄)
end

function cheb_xrange(x̃range,x̃)
        imax = _searchsortedfirst(x̃range, x̃)
        imin = imax - 1
        x̃minᵢ = x̃range[imin]
        x̃maxᵢ = x̃range[imax]
        i = imin
    x̄ = (2*x̃ - (x̃maxᵢ + x̃minᵢ)) / (x̃maxᵢ - x̃minᵢ)
    return x̄,i
end

function cheb_eval(cheb::AbstractVector{T},x::S) where {T,S}
    R = promote_type(T, S)
    l = length(cheb)
    l == 0 && return zero(R)
    l == 1 && return R(cheb[1])
    c0 = cheb[l - 1]
    c1 = cheb[l]
    for i in (l-2):-1:1
        c0, c1 = cheb[i] - c1, c0 + c1 * 2x
    end
    return R(c0 + c1 * x)
end
