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
saturation_pressure(model,Tsat,SuperAncSaturation()) #solve using cubic superancillary
```

Bubble point pressure
```julia
model = PCSAFT(["methanol","cyclohexane"])
T = 313.15
z = [0.5,0.5]
bubble_pressure(model,T,z) #using default method (chemical potential equality)
bubble_pressure(model,T,z,FugBubblePressure(y0 = = [0.6,0.4], p0 = 5e4)) #using isofugacity criteria with starting points
```
"""
abstract type ThermodynamicMethod end
Base.show(io::IO,::MIME"text/plain",x::ThermodynamicMethod) = show_as_namedtuple(io,x)
Base.show(io::IO,x::ThermodynamicMethod) = show_as_namedtuple(io,x)

function NLSolvers.NEqOptions(method::ThermodynamicMethod)
    return NEqOptions(f_limit = method.f_limit,
                    f_abstol = method.atol,
                    f_reltol = method.rtol,
                    maxiter = method.max_iters)
end

mw(model::EoSModel) = model.params.Mw.values

group_Mw(model::EoSModel) = group_Mw(model.params.Mw.values,model.groups)
function group_Mw(Mw_gc::SingleParam,groups::GroupParam)
    n = length(groups.components)
    mw_comp = zeros(eltype(Mw_gc.values),n)
    v = groups.n_flattenedgroups
    mw_gc = Mw_gc.values
    for i in 1:n
        mw_comp[i] = dot(mw_gc,v[i])
    end
    return mw_comp
end

function group_molecular_weight(groups::GroupParameter,mw,z = @SVector [1.])
    res = zero(first(z))
    for i in 1:length(groups.n_flattenedgroups)
        ni = groups.n_flattenedgroups[i]
        mwi = dot(ni,mw)
        res +=z[i]*mwi
    end
    return 0.001*res/sum(z)
end

comp_molecular_weight(mw,z = SA[1.0]) = 0.001*dot(mw,z)
molecular_weight(model) = molecular_weight(model,SA[1.0])
molecular_weight(model::EoSModel,z) = __molecular_weight(model,z)
molecular_weight(mw::AbstractVector,z) = comp_molecular_weight(mw,z)
molecular_weight(mw::SingleParam,z) = comp_molecular_weight(mw.values,z)

function __molecular_weight(model,z)
    MW = mw(model)
    if has_groups(model)
        return group_molecular_weight(model.groups,MW,z)
    else
        return comp_molecular_weight(MW,z)
    end
end

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
is_supercritical(str::String) = is_supercritical(Symbol(str))

const UNKOWN_STR = (:unknown,:UNKNOWN)

"""
    is_unknown(x::Union{Symbol,String})

Returns `true` if the symbol is in `(:unknown,:UNKNOWN)`.

If a string is passed, it is converted to symbol.
"""
is_unknown(sym::Symbol) = sym in UNKOWN_STR
is_unknown(str::String) = is_unknown(Symbol(str))

const SOLID_STR = (:solid,:SOLID,:s)
"""
    is_solid(x::Union{Symbol,String})
    is_solid(x::EoSModel)

Returns `true` if the symbol is in `(:solid,:SOLID,:s)`.
if `x` is an `EoSModel`, it will return if the model is able to contain a solid phase. In this case, defaults to `false`
If a string is passed, it is converted to symbol.
"""
is_solid(sym::Symbol) = sym in SOLID_STR
is_solid(str::String) = is_solid(Symbol(str))
is_solid(model::EoSModel) = false


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
    hastype = (Base.IteratorEltype(typeof(iterator)) === Base.HasEltype()) && (eltype(iterator) !== Any)
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

∑(x::AbstractArray) = sum(x)
∑(f,x::AbstractArray) = sum(f,x)

function ∑(fn,iterator)
    len = Base.IteratorSize(typeof(iterator)) === Base.HasLength()
    hastype = (Base.IteratorEltype(typeof(iterator)) === Base.HasEltype()) && (eltype(iterator) !== Any)
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

function is_ad_input(model,V,T,z)
    #model_primal = Solvers.primal_eltype(model)
    #TODO: do something if the model parameters themselves use ForwardDiff 
    V_primal,T_primal,z_primal = Solvers.primalval(V),Solvers.primalval(T),Solvers.primalval(z)
    type = Base.promote_eltype(V,T,z)
    primal_type = Base.promote_eltype(V_primal,T_primal,z_primal)
    return primal_type != type
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
    μ = Base.promote_eltype(V,T,z)
    return StaticArrays.similar_type(z,μ)
end

function gradient_type(V,T,z::AbstractVector)
    μ = Base.promote_eltype(V,T,z)
    return Vector{μ}
end

function gradient_type(V,T,z::FractionVector{TT,UU}) where {TT,UU<:AbstractVector}
    μ = Base.promote_eltype(V,T,z)
    return Vector{μ}
end

#= todo: fix this
function gradient_type(V,T,z::FractionVector{TT,TT}) where {TT}
    μ = Base.promote_eltype(V,T,z)
    return SVector{2, μ}
end =#



"""
    init_preferred_method(method,model,kwargs)

Returns the preferred method for a combination of model and function, with the specified kwargs.

"""
function init_preferred_method(method,model) end

"""
    get_k(model)::VarArg{Matrix}

Returns a matrix of "k-values" binary interaction parameters used by the input `model`. Returns `nothing` if the model cannot return the k-values matrix.
In the case of multiple k-values (as is the case in T-dependent values, i.e: k(T) = k1 + k2*T), it will return a tuple of matrices corresponding to each term in the k-value expression.
Note that some models do not store the k-value matrix directly, but they contain the value in an indirect manner. for example, cubic EoS store `a[i,j] = f(a[i],a[j],k[i,j])`, where `f` depends on the mixing rule.
"""
get_k(model::EoSModel) = nothing

"""
    get_l(model)::VarArg{Matrix}

returns a matrix of "l-values" binary interaction parameters used by the input `model`. Returns `nothing` if the model cannot return the l-values matrix.
In the case of multiple l-values (as is the case in T-dependent values, i.e: l(T) = l1 + l2*T), it will return a tuple of matrices corresponding to each term in the l-value expression.
Note that some models do not store the l-value matrix directly, but they contain the value in an indirect manner. for example, cubic EoS store `b[i,j] = f(b[i],b[j],l[i,j])`, where `f` depends on the mixing rule.
"""
get_l(model::EoSModel) = nothing

"""
    set_k!(model,k)
    set_k!(model,ki...)

Sets the model "k-values" binary interaction parameter to the input matrix `k`. If the input model requires multiple k-matrices (as is the case for T-dependent values, i.e: k(T) = k1 + k2*T), then you must call `set_k!` with all the matrices as input (`set_k!(model,k1,k2)`).

"""
set_k!(model::EoSModel,k) = throw(ArgumentError("$(typeof(model)) does not have support for setting k-values"))

"""
    set_l!(model,l)
    set_l!(model,li...)

Sets the model "l-values" binary interaction parameter to the input matrix `l`. If the input model requires multiple l-matrices (as is the case for T-dependent values, i.e: l(T) = l1 + l2*T), then you must call `set_l!` with all the matrices as input (`set_l!(model,l1,l2)`).

"""
set_l!(model::EoSModel,k) = throw(ArgumentError("$(typeof(model)) does not have support for setting l-values"))

export get_k,set_k!
export get_l,set_l!

include("initial_guess.jl")
include("differentials.jl")
include("VT.jl")
include("isochoric.jl")
include("fugacity_coefficient.jl")
include("property_solvers/property_solvers.jl")
include("tpd.jl")
include("stability.jl")
include("pT.jl")
include("property_solvers/Tproperty.jl")
include("property_solvers/Pproperty.jl")
include("XY_methods/VT.jl")
include("XY_methods/PS.jl")
include("XY_methods/TS.jl")
include("XY_methods/PH.jl")
include("XY_methods/QX.jl")

include("property_solvers/spinodal.jl")
