abstract type CPAAlphaModel <: AlphaModel end

struct CPAAlphaParam <: EoSParam
    c1::SingleParam{Float64}
end

@newmodelsimple CPAAlpha CPAAlphaModel CPAAlphaParam
export CPAAlpha

"""
    CPAAlpha <: CPAAlphaModel
    
    CPAAlpha(components;
    userlocations=String[],
    verbose::Bool=false)

## Input Parameters

- `c1`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for `CPA` EoS.
```
αᵢ = (1+c¹ᵢ(1-√(Trᵢ)))^2

```

"""
CPAAlpha
default_locations(::Type{CPAAlpha}) = ["SAFT/CPA/CPA_like.csv"]

function α_function(model::CubicModel,V,T,z,alpha_model::CPAAlphaModel)
    Tc = model.params.Tc.values
    c1  = alpha_model.params.c1.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        Tr = T/Tc[i]
        α[i] = (1+c1[i]*(1-√(Tr)))^2
    end
    return α
end
