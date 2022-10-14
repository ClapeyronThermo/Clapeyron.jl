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

function saturation_pressure(model::SaturationModel,T,method::SaturationMethod)
    single_component_check(saturation_pressure,model)
    T = T*(T/T)
    return saturation_pressure_impl(model,T,SaturationCorrelation())
end

function saturation_temperature(model::SaturationModel,p,method::SaturationMethod)
    single_component_check(saturation_temperature,model)
    p = p*(p/p)
    return saturation_temperature_impl(model,p,SaturationCorrelation())
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
    f0(T) = first(saturation_pressure_impl(model,T,method)) - p
    prob = Roots.ZeroProblem(f0,T0)
    sol = Roots.solve(prob)
    return sol,nan,nan
end

#=
    tc = temperature(model,CriticalPoint())
    pc = pressure(model,CriticalPoint())
    t7 = 0.7*tc
    p7 = pressure_impl(QuickStates.sat_t(),model,t7)
    
    h = 2.3333333333333335*log(pc/p7)
    return 1/(1-log(p/pc)/h)*tc
end
=#
eos(model,V,T,z=SA[1.0]) = not_eos_error(model)

include("LeeKeslerSat/LeeKeslerSat.jl")
include("DIPPR101Sat/DIPPR101Sat.jl")
