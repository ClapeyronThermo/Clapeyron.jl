# extract all sites
const LIQUID_STR = (:liquid,:LIQUID,:L,:l)


"""
    is_liquid(x::Union{Symbol,String})   

Returns `true` if the symbol is in `(:liquid,:LIQUID,:L,:l)`.

If a string is passed, it is converted to symbol.
"""
is_liquid(sym::Symbol) = sym in LIQUID_STR
is_liquid(str::String) = is_liquid(Symbol(str)) 

const VAPOUR_STR = (:vapor,:VAPOR,:VAPOUR,:vapour,:g,:G,:v,:V,:gas,:GAS)

"""
    is_vapour(x::Union{Symbol,String})   

Returns `true` if the symbol is in `(:vapor,:VAPOR,:VAPOUR,:vapour,:g,:G,:v,:V,:gas,:GAS)`.

If a string is passed, it is converted to symbol.
"""
is_vapour(sym::Symbol) = sym in VAPOUR_STR
is_vapour(str::String) = is_vapour(Symbol(str))

const SUPERCRITICAL_STR = (:sc,:SC,:supercritical,:SUPERCRITICAL)


"""
    is_supercritical(x::Union{Symbol,String})   

Returns `true` if the symbol is in `(:sc,:SC,:supercritical,:SUPERCRITICAL)`.

If a string is passed, it is converted to symbol.
"""
is_supercritical(sym::Symbol) = sym in SUPERCRITICAL_STR
is_supercritical(str::String) = is_vapour(Symbol(str))


"""
    ∑(iterator)

equivalent to `sum(iterator,init=0.0)`. 

"""
function ∑(iterator) #not collecting is faster
    return reduce(Base.add_sum,iterator,init=0.0)
end

struct FractionVector{T,V} <: AbstractVector{T}
    vec::V
    val::T
end
#=

Fraction Vector
useful when expressing fractions in n-1 components.
the last component is calculated at build time.
it allocates less than creating a new vector or appending.
=#
##
function FractionVector(v::AbstractVector)
    a = one(eltype(v))
    for vi in v
        vi < 0 && throw(DomainError(vi,"all elements of a fraction vector should be positive."))
        a -= vi
    end
    a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
    return FractionVector(v,a) 
end

function FractionVector(v::Real)
    return FractionVector(SA[v])
end


@inline Base.eltype(v::FractionVector{T}) where T = T
@inline Base.length(v::FractionVector)::Int = Base.length(v.vec) + 1
@inline Base.size(v::FractionVector) = (length(v),)
@inline function Base.getindex(v::FractionVector,i)
    if length(v) == i
        return v.val
    else
        return v.vec[i]
    end
end

Base.IndexStyle(::Type{<:FractionVector}) = IndexLinear()