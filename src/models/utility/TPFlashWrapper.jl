#wrapper used to cache results in case of activity models and CompositeModel
struct PTFlashWrapper{T,S} <: EoSModel
    components::Vector{String}
    model::T
    sat::Vector{S}
    fug::Vector{S}
    μ::Vector{S}
end

Base.length(model::PTFlashWrapper) = length(model.model)

function tp_flash_K0(wrapper::PTFlashWrapper,p,T)
    first.(wrapper.sat) ./ p
end

function tp_flash_K0!(K,wrapper::PTFlashWrapper,p,T)
    K .=  first.(wrapper.sat) ./ p 
end

function PTFlashWrapper(model::EoSModel,p,T::Number,equilibrium::Symbol)
    pures = split_model(model)
    RT = R̄*T
    sats = saturation_pressure.(pures,T)
    vv_pure = last.(sats)
    p_pure = first.(sats)
    μpure = only.(VT_chemical_potential_res.(__gas_model.(pures),vv_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    g_pure = [VT_gibbs_free_energy(__gas_model(pures[i]),vv_pure[i],T) for i in 1:length(model)]
    return PTFlashWrapper(model.components,model,sats,ϕpure,g_pure)
end
