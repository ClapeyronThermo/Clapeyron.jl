abstract type SaturationModel <: EoSModel end

struct SaturationCorrelation <: SaturationMethod end

function saturation_pressure(model::SaturationModel,T,method::SaturationMethod)
    single_component_check(saturation_pressure,model)
    T = T*(T/T)
    return saturation_pressure_impl(model,T,SaturationCorrelation())
end

eos(model,V,T,z=SA[1.0]) = not_eos_error(model)

include("LeeKeslerSat/LeeKeslerSat.jl")
