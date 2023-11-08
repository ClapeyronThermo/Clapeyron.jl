#wrapper used to cache results in case of activity models and CompositeModel
struct PTFlashWrapper{T,S,F,M} <: EoSModel
    components::Vector{String}
    model::T
    sat::Vector{S}
    fug::Vector{F}
    μ::Vector{M}
end

Base.length(model::PTFlashWrapper) = length(model.model)

function tp_flash_K0(wrapper::PTFlashWrapper,p,T)
    first.(wrapper.sat) ./ p
end

function tp_flash_K0!(K,wrapper::PTFlashWrapper,p,T)
    K .=  first.(wrapper.sat) ./ p 
end
