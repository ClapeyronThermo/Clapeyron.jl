mw(model::EoSModel) = model.params.Mw.values

function group_molecular_weight(groups::GroupParam,mw,z = @SVector [1.])
    n = groups.n_flattenedgroups
    res = zero(first(z))
    Σz = sum(z)
    @inbounds for i in 1:length(groups.components)
        ni = n[i]
        gi = groups.i_groups[i]
        mwi = zero(res)
        for idx in 1:length(gi)
            j = gi[idx]
            mwi += mw[j]*ni[j]
        end
        res +=z[i]*mwi
    end
    return 0.001*res/Σz
end

comp_molecular_weight(mw,z = @SVector [1.]) = 0.001*dot(mw,z)

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
function ∑(iterator)
    len = Base.IteratorSize(typeof(iterator)) === Base.HasLength()
    hastype =  (Base.IteratorEltype(typeof(iterator)) === Base.HasEltype()) && (eltype(iterator) !== Any)
    local _0
    if hastype
        _0 = zero(eltype(iterator))
    else
        _0 = 0.0
    end
    len && iszero(length(iterator)) && return _0
    !len && return reduce(Base.add_sum,iterator,init=_0)
    return sum(iterator)
end

function ∑(fn,iterator) 
    len = Base.IteratorSize(typeof(iterator)) === Base.HasLength()
    hastype =  (Base.IteratorEltype(typeof(iterator)) === Base.HasEltype()) && (eltype(iterator) !== Any)
    local _0
    if hastype
        _0 = zero(eltype(iterator))
    else
        _0 = 0.0
    end
    len && iszero(length(iterator)) && return _0
    !len && return mapreduce(fn,Base.add_sum,iterator,init=_0)
    return sum(fn,iterator)
end

"""
    xlogx(x::Real)
Return `x * log(x)` for `x ≥ 0`, handling ``x = 0`` by taking the downward limit.

copied from LogExpFunctions.jl
"""
function xlogx(x::Real)
    _0 = zero(x)
    iszero(x) && return _0
    ifelse(x >= _0,x*Base.log(max(_0,x)),_0/_0)
end

@inline function nan_num(V,T,z)
    _0 = zero(V+T+first(z))
    _0/_0
end

@inline function negative_vt(V,T)::Bool
    _0 = zero(V+T)
    (T <= _0) | (V <= _0)   
end

#the only macro needed in methods
"""
    @nan(function_call,default=NaN)

Wraps the function in a `try-catch` block, and if a `DomainError` or `DivideError` is raised, then returns `default`.
for better results, its best to generate the default result beforehand
"""
macro nan(Base.@nospecialize(fcall),default = nothing)
    quote
      try $fcall
      catch err
        if err isa Union{DomainError,DivideError}
          $default
        else
          rethrow(err)
        end
      end
    end |> esc
end

function gradient_type(V,T,z::StaticArray)
    μ = typeof(V+T+first(z))
    return StaticArrays.similar_type(z,μ)
end

function gradient_type(V,T,z::Vector)
    μ = typeof(V+T+first(z))
    return Vector{μ}
end

function gradient_type(V,T,z::FractionVector)
    μ = typeof(V+T+first(z))
    return Vector{μ}
end

include("initial_guess.jl")
include("differentials.jl")
include("VT.jl")
include("property_solvers/property_solvers.jl")
include("pT.jl")
include("unitful_base.jl")
include("unitful_methods.jl")

