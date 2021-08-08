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
function ∑ end

if Base.VERSION < v"1.6"
    function ∑(iterator)
        _0 = 0.0
        return reduce(Base.add_sum,iterator,init=_0)
    end
    ∑(iterator) = sum(iterator)
    function ∑(iterator::AbstractArray) 
        if iszero(length(iterator))
            return zero(eltype(iterator))
        else
            return sum(iterator)
        end
    end

    function ∑(fn,iterator) #not collecting is faster
        #_0 = ifelse(length(iterator)>0,zero(first(iterator)),0.0)
        _0 = 0.0
        return mapreduce(fn,Base.add_sum,iterator,init=_0)
    end
else
    ∑(iterator) = sum(iterator)
    function ∑(iterator::AbstractArray) 
        ifelse(iszero(length(iterator)),zero(eltype(iterator)),sum(iterator))
    end
    ∑(fn,iterator) = sum(fn,iterator) 
end




"""
    xlogx(x::Number)
Return `x * log(x)` for `x ≥ 0`, handling ``x = 0`` by taking the downward limit.

copied from LogExpFunctions.jl
"""
function xlogx(x::Number)
    result = x * log(x)
    ifelse(iszero(x), zero(result), result)
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
    any(x->x<0,v) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
    a -=sum(v)
    a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
    return FractionVector(v,a) 
end

function FractionVector(v::Real)
    a = one(v) 
    (v < zero(v)) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
    a -= v
    a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
    return FractionVector(v,a)
end


@inline Base.eltype(v::FractionVector{T}) where T = T
@inline Base.length(v::FractionVector)::Int = Base.length(v.vec) + 1

@inline function Base.length(v::FractionVector{T,<:Real})::Int where T
    return 2
end

@inline Base.size(v::FractionVector) = (length(v),)
@inline function Base.getindex(v::FractionVector,i::Int)
    @boundscheck checkbounds(v, i)
    return ifelse(length(v)==i,v.val,v.vec[min(length(v.vec),i)])
end

@inline function Base.getindex(v::FractionVector{T,<:Real},i::Int) where T
    @boundscheck checkbounds(v, i)
    return ifelse(i==1,v.vec,v.val)
end

Base.IndexStyle(::Type{<:FractionVector}) = IndexLinear()

unicode_subscript(a::Int) = join('₀'+d for d in reverse(digits(a)))

"""
    split_model(model::EoSModel)
Takes in a model for a multi-component system and returns a vector of model for each pure system.
"""
function split_model(model::EoSModel)
    pure = []
    for i ∈ model.components
        append!(pure,[eval(nameof(typeof(model)))([i])])
    end
    return pure
end