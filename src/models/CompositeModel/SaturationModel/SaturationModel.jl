"""
    SaturationModel <: EoSModel

Abstract type for Saturation correlation models.
"""
abstract type SaturationModel <: EoSModel end

"""
    SaturationCorrelation <: SaturationMethod

Saturation method used for dispatch on saturation correlations.
"""
struct SaturationCorrelation <: SaturationMethod end

function init_preferred_method(method::typeof(saturation_pressure),model::SaturationModel,kwargs)
    return SaturationCorrelation()
end

function init_preferred_method(method::typeof(saturation_temperature),model::SaturationModel,kwargs)
    return SaturationCorrelation()
end
 
function saturation_temperature_impl(model::SaturationModel,p,method::SaturationCorrelation)
    if any(has_dual,(model,p))
        λmodel,λp = primalval(model),primalval(p)
        λT0 = saturation_temperature_corr_x0(λmodel,λp)
        λf0(T) = log(first(saturation_pressure_impl(λmodel,T,method))/λp)
        λprob = Roots.ZeroProblem(λf0,λT0)
        λT = Roots.solve(λprob)
        sol = saturation_temperature_corr_ad(λT,tup,tup_primal)
    else
        T0 = saturation_temperature_corr_x0(model,p)
        f0(T) = log(first(saturation_pressure_impl(model,T,method))/p)
        prob = Roots.ZeroProblem(f0,T0)
        sol = Roots.solve(prob)
    end
    nan = zero(sol)/zero(sol)
    return sol,nan,nan
end

function saturation_temperature_corr_x0(model,p)
    Tc,Pc,_ = crit_pure(model)
    _1 = one(Base.promote_eltype(model,p))
    nan = zero(_1)/zero(_1)
    p > Pc && return nan
    #critical interpolation
    T_07 = 0.7*Tc
    p_07,_,_ = saturation_pressure_impl(model,T_07,SaturationCorrelation())
    h = 2.3333333333333335*log(Pc/p_07)
    T0 = _1/(1-log(p/Pc)/h)*Tc
    return T0
end

function saturation_temperature_corr_ad(λT,tup,tup_primal)
    f(T,tups) = begin
        model,p = tups
        return first(saturation_pressure_impl(model,T,SaturationCorrelation())) - p
    end
    ∂T,∂vl,∂vv = __gradients_for_root_finders(λT,tup,tup_primal,f)
    return ∂T,∂vl,∂vv
end

function saturation_temperature_ad(__model::SaturationModel,result,tup,tup_primal)
    ∂T = saturation_temperature_corr_ad(first(result),tup,tup_primal)
    nan = zero(∂T)/zero(∂T)
    return ∂T,nan,nan
end

#this method allows to use a Saturation Model as a whole fluid model. it supposes ideal gas and no info about the liquid phase (ZeroLiquid)
function init_puremodel(model::SaturationModel,components,userlocations,verbose)
    _components = format_components(components)
    fluid = CompositeModel(_components,gas=BasicIdeal(),liquid=ZeroLiquid(),saturation = model)
    return EoSVectorParam(fluid,_components)
end 

function dpdT_saturation(model::SaturationModel,v1::Number,v2,T)
    sat(_T) = first(saturation_pressure(model,_T))
    return Solvers.derivative(sat,T)
end

include("LeeKeslerSat/LeeKeslerSat.jl")
include("DIPPR101Sat/DIPPR101Sat.jl")
include("AntoineSat/AntoineEqSat.jl")
include("PolExpSat/PolExpSat.jl")
