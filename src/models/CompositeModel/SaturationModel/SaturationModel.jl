"""
    SaturationModel <: EoSModel

Abstract type for Saturation correlation models.
"""
abstract type SaturationModel <: EoSModel end

"""
    SaturationCorrelation <: SaturationMethod

saturation method used for dispatch on saturation correlations.
"""
struct SaturationCorrelation <: SaturationMethod end

function init_preferred_method(method::typeof(saturation_pressure),model::SaturationModel,kwargs)
    return SaturationCorrelation()
end

function init_preferred_method(method::typeof(saturation_temperature),model::SaturationModel,kwargs)
    return SaturationCorrelation()
end
 
function saturation_temperature_impl(model::SaturationModel,p,method::SaturationCorrelation)
    Tc,Pc,_ = crit_pure(model)
    nan = zero(p)/zero(p)
    p > Pc && return (nan,nan,nan)
    #critical interpolation
    T_07 = 0.7*Tc
    p_07,_,_ = saturation_pressure_impl(model,T_07,method)
    h = 2.3333333333333335*log(Pc/p_07)
    T0 = 1/(1-log(p/Pc)/h)*Tc
    f0(T) = log(first(saturation_pressure_impl(model,T,method))/p)
    prob = Roots.ZeroProblem(f0,T0)
    sol = Roots.solve(prob)
    return sol,nan,nan
end

eos(model,V,T,z=SA[1.0]) = not_eos_error(model)

include("LeeKeslerSat/LeeKeslerSat.jl")
include("DIPPR101Sat/DIPPR101Sat.jl")
include("PolExpSat/PolExpSat.jl")
