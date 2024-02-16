abstract type EmpiricIdealModel <: IdealModel end

struct EmpiricIdeal <: EmpiricIdealModel
    components::Vector{String}
    params::MultiFluidParam
    pures::Vector{SingleFluidIdeal}
    Rgas::Float64
    references::Vector{String}
end

struct EmpiricIdealfromMulti{𝔸} <: EmpiricIdealModel
    components::Vector{String}
    params::MultiFluidParam
    pures::Vector{SingleFluid{𝔸}}
    Rgas::Float64
    references::Vector{String}
end

Rgas(m::EmpiricIdealModel) = m.Rgas

"""
    EmpiricIdeal(components;
    pure_userlocations = String[],
    estimate_pure = false,
    coolprop_userlocations = true,
    Rgas = R̄,
    verbose = false)

## Input parameters
- JSON data (CoolProp and teqp format)

## Description

Instantiates the ideal part of a multi-component Empiric EoS model. `Rgas` can be used to set the value of the gas constant that is used during property calculations.

If `coolprop_userlocations` is true, then Clapeyron will try to look if the fluid is present in the CoolProp library.

If `estimate_pure` is true, then, if a JSON is not found, the pure model will be estimated, using the `XiangDeiters` model

"""
EmpiricIdeal

function EmpiricIdeal(components;
    userlocations = String[],
    coolprop_userlocations = true,
    Rgas = R̄,
    reference_state = nothing,
    verbose = false)

    components = format_components(components)
    pures = [
        SingleFluidIdeal(comp;
        userlocations = userlocations,
        verbose = verbose,
        coolprop_userlocations = coolprop_userlocations,
        Rgas = Rgas
        )
        for comp in components]
    params = MultiFluidParam(components,pures,reference_state)
    references = unique!(reduce(vcat,pure.references for pure in pures))
    model = EmpiricIdeal(components,params,pures,Rgas,references)
    set_reference_state!(model,verbose = verbose)
    return model
end

function idealmodel(m::MultiFluid)
    EmpiricIdealfromMulti(m.components,m.params,m.pures,m.Rgas,m.references)
end

set_reference_state!(model::EmpiricIdealModel;verbose = false) = set_reference_state_empiric!(model;verbose)

function a_ideal(model::EmpiricIdealModel,V,T,z,∑z = sum(z))
    #log(δi) = log(ρ * vc[i]) = -log(V) + log(sum(z)) + log(vc[i])
    res = zero(V+T+first(z))
    m₀ = model.pures
    Tinv = 1/T
    Tc = model.params.Tc
    vc = model.params.Vc
    for i in 1:length(model)
        m₀ᵢ = m₀[i]
        a₀ᵢ = __get_k_alpha0(m₀ᵢ)*reduced_a_ideal(m₀ᵢ,Tc[i] * Tinv)
        zᵢ = z[i]
        res += zᵢ*a₀ᵢ
        res += xlogx(zᵢ,vc[i])
        res
    end
    res /= ∑z
    res -= log(V)
    return res
end

export EmpiricIdeal


