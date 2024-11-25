abstract type RKAlphaModel <: AlphaModel end

@newmodelsingleton RKAlpha RKAlphaModel

"""
    RKAlpha <: RKAlphaModel
    
    RKAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- none

## Description

Cubic alpha `(α(T))` model. Default for [`RK`](@ref) EoS.
```
αᵢ = 1/√(Trᵢ)
Trᵢ = T/Tcᵢ
```

## Model Construction Examples
```
# Because this model does not have parameters, all those constructors are equivalent:
alpha = RKAlpha()
alpha = RKAlpha("water")
alpha = RKAlpha(["water","carbon dioxide"])
```
"""
RKAlpha

function α_function!(α,model::CubicModel,alpha_model::RKAlphaModel,T)
    Tc = model.params.Tc.values
    for i in @comps
        Tr = T/Tc[i]
        α[i] = 1/√(Tr)
    end
    return α
end

function α_function(model::CubicModel,alpha_model::RKAlphaModel,T,i::Int)
    Tc = model.params.Tc.values[i]
    Tr = T/Tc
    α = 1/√(Tr)
end