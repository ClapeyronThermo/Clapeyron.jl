"""
    GammaPhi{γ,Φ} <: EoSModel

wrapper struct to signal that a `CompositeModel` uses an activity model in conjunction with a fluid.
"""
struct GammaPhi{γ,Φ} <: EoSModel
    components::Vector{String}
    activity::γ
    fluid::Φ
end

function init_preferred_method(method::typeof(saturation_pressure),model::GammaPhi,kwargs)
    return init_preferred_method(saturation_pressure,model.fluid,kwargs)
end

function init_preferred_method(method::typeof(saturation_temperature),model::GammaPhi,kwargs)
    return init_preferred_method(saturation_temperature,model.fluid,kwargs)
end

function saturation_pressure(model::GammaPhi,T,method::SaturationMethod)
    return saturation_pressure(model.fluid,T,method)
end

crit_pure(model::GammaPhi) = crit_pure(model.fluid)
x0_sat_pure(model::GammaPhi,T) = x0_sat_pure(model.fluid)
x0_psat(model::GammaPhi,T) = x0_psat(model.fluid,T)

function saturation_temperature(model::GammaPhi,p,method::SaturationMethod)
    return saturation_temperature(model.fluid,p,method)
end

function init_preferred_method(method::typeof(tp_flash),model::GammaPhi,kwargs)
    RRTPFlash(;kwargs...)
end