#wrapper used to cache results in case of activity models and CompositeModel
struct PTFlashWrapper{T,R,S} <: EoSModel
    components::Vector{String}
    model::T
    sat::Vector{R}
    fug::Vector{S}
    μ::Vector{S}
    equilibrium::Symbol
end

Base.length(model::PTFlashWrapper) = length(model.model)

function tp_flash_K0(wrapper::PTFlashWrapper,p,T,z)
    first.(wrapper.sat) ./ p
end

function tp_flash_K0!(K,wrapper::PTFlashWrapper,p,T,z)
    K .=  first.(wrapper.sat) ./ p 
end

function PTFlashWrapper(model::EoSModel,p,T::Number,equilibrium::Symbol)
    pures = split_pure_model(model)
    RT = R̄*T
    sats = saturation_pressure.(pures,T)
    vv_pure = last.(sats)
    p_pure = first.(sats)
    μpure = only.(VT_chemical_potential_res.(gas_model.(pures),vv_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    g_pure = [VT_gibbs_free_energy(gas_model(pures[i]),vv_pure[i],T) for i in 1:length(model)]
    return PTFlashWrapper(component_list(model),model,sats,ϕpure,g_pure,equilibrium)
end

gas_model(model::PTFlashWrapper) = gas_model(model.model)

function volume_impl(model::PTFlashWrapper, p, T, z, phase, threaded, vol0)
    volume_impl(model.model, p, T, z, phase, threaded, vol0)
end