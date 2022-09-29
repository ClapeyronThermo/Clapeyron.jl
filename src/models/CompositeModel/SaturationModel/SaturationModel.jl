abstract type SaturationModel <: EoSModel end

struct SaturationCorrelation <: SaturationMethod end

function saturation_pressure(model::SaturationModel,T)
    saturation_pressure(model,T,SaturationCorrelation())
end

function saturation_pressure_impl(model::SaturationModel,T,method::SaturationMethod)
    throw(error("$method not supported by Saturation correlation models"))
end

include("LeeKeslerSat/LeeKeslerSat.jl")
