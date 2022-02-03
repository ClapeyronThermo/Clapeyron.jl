abstract type SoaveAlphaModel <: AlphaModel end

struct SoaveAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple SoaveAlpha SoaveAlphaModel SoaveAlphaParam
export SoaveAlpha

"""
    SoaveAlpha <: SoaveAlphaModel
    
    SoaveAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `w`: Single Parameter (`Float64`)

## Model Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for `SRK` EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.480 + 1.547ωᵢ - 0.176ωᵢ^2
```
to use different polynomial coefficients for `mᵢ`, overload `Clapeyron.α_m(::CubicModel,::SoaveAlphaModel) = (c₁,c₂,...cₙ)`

"""
SoaveAlpha

function SoaveAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = SoaveAlphaParam(acentricfactor)
    model = SoaveAlpha(packagedparams, verbose=verbose)
    return model
end

@inline α_m(model,::SoaveAlpha) = (0.480,1.547,-0.176)

function α_function(model::CubicModel,V,T,z,alpha_model::SoaveAlphaModel)
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    α = zeros(typeof(T),length(Tc))
    coeff = α_m(model,alpha_model)
    for i in @comps
        ωi = ω[i]
        Tr = T/Tc[i]
        m = evalpoly(ωi,coeff)
        α[i] = (1+m*(1-√(Tr)))^2
    end
    return α
end