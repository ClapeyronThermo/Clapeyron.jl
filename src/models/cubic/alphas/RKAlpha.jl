abstract type RKAlphaModel <: AlphaModel end

@newmodelsingleton RKAlpha RKAlphaModel

"""
    RKAlpha <: RKAlphaModel
    
    RKAlpha(components;
    userlocations=String[],
    verbose::Bool=false)

## Input Parameters

- `w`: Single Parameter (`Float64`)

## Model Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`RK`](@ref) EoS.
```
αᵢ = 1/√(Trᵢ)
Trᵢ = T/Tcᵢ
```

"""
RKAlpha

function α_function(model::CubicModel,V,T,z,alpha_model::RKAlphaModel)
    Tc = model.params.Tc.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        Tr = T/Tc[i]
        α[i] = 1 /√(Tr)
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::RKAlphaModel)
    Tc = model.params.Tc.values[1]
    Tr = T/Tc
    α = 1 /√(Tr)
end
