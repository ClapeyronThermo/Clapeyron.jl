"""
    ThermodynamicMethod

Abstract type for all thermodynamic methods.

normally, a thermodynamic method has the form: `property(model,state..,method::ThermodynamicMethod)`.
All methods used in this way subtype `ThermodynamicMethod`.

## Examples
Saturation pressure:
```julia
model = PR(["water"])
Tsat = 373.15
saturation_pressure(model,Tsat) #using default method (chemical potential with volume base)
saturation_pressure(model,Tsat,SuperAncSaturation()) #solve using cubic superancilliary
```

Bubble point pressure
```julia
model = PCSAFT(["methanol","cyclohexane"])
T = 313.15
z = [0.5,0.5]
bubble_pressure(model,T,z) #using default method (chemical potential equality)
bubble_pressure(model,T,z,FugBubblePressure(y0 =  = [0.6,0.4], p0 = 5e4)) #using isofugacity criteria with starting points
```
"""
abstract type ThermodynamicMethod end

function NLSolvers.NEqOptions(method::ThermodynamicMethod)
    return NEqOptions(f_limit = method.f_limit,
                    f_abstol = method.atol,
                    f_reltol = method.rtol,
                    maxiter = method.max_iters)
end

mw(model::EoSModel) = model.params.Mw.values

function group_molecular_weight(groups::GroupParam,mw,z = @SVector [1.])
    res = zero(first(z))
    for ni in groups.n_flattenedgroups
        mwi = dot(ni,mw)
        res +=z[i]*mwi
    end
    return 0.001*res/sum(z)
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

const SOLID_STR = (:solid,:SOLID,:s)
"""
    is_solid(x::Union{Symbol,String})

Returns `true` if the symbol is in `(:solid,:SOLID,:s)`.

If a string is passed, it is converted to symbol.
"""
is_solid(sym::Symbol) = sym in SOLID_STR
is_solid(str::String) = is_vapour(Symbol(str))


const VLE_STR = (:vle,:lve,:vl,:lv)
"""
    is_vle(x::Union{Symbol,String})

Returns `true` if the symbol is in `(:vle,:lve,:vl,:lv)`.

If a string is passed, it is converted to symbol.
"""
is_vle(sym::Symbol) = sym in VLE_STR
is_vle(str::String) = is_vle(Symbol(str))

const LLE_STR = (:lle,:ll)
"""
    is_lle(x::Union{Symbol,String})

Returns `true` if the symbol is in `(:lle,:ll)`.

If a string is passed, it is converted to symbol.
"""
is_lle(sym::Symbol) = sym in LLE_STR
is_lle(str::String) = is_lle(Symbol(str))

function canonical_phase(phase::Symbol)
     if is_liquid(phase)
        return :liquid
     elseif is_vapour(phase)
        return :vapour
     else
        return phase
     end
end

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
include("fugacity_coefficient.jl")
include("property_solvers/property_solvers.jl")
include("tpd.jl")
include("stability.jl")
include("pT.jl")
include("unitful_base.jl")
include("unitful_methods.jl")
