abstract type ClausiusAlphaModel <: AlphaModel end
@newmodelsingleton ClausiusAlpha ClausiusAlphaModel
export ClausiusAlpha

"""
    ClausiusAlpha <: ClausiusAlphaModel
    
    ClausiusAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- none

## Description

Cubic alpha `(α(T))` model. Default for [`Clausius`](@ref)  and [`Berthelot`]
```
αᵢ = 1/Trᵢ
Trᵢ = T/Tcᵢ
```

## Model Construction Examples
```
# Because this model does not have parameters, all those constructors are equivalent:
alpha = ClausiusAlpha()
alpha = ClausiusAlpha("water")
alpha = ClausiusAlpha(["water","carbon dioxide"])
```
"""
ClausiusAlpha

function α_function(model::CubicModel,V,T,z,alpha_model::ClausiusAlphaModel)
    Tc = model.params.Tc.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        α[i] = Tc[i]/T
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::ClausiusAlphaModel)
    Tc = model.params.Tc.values[1]
    α = Tc/T
end