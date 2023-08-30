abstract type sCPAAlphaModel <: CPAAlphaModel end


@newmodelsimple sCPAAlpha sCPAAlphaModel CPAAlphaParam

"""
    sCPAAlpha <: sCPAAlphaModel
    
    CPAAlpha(components;
    userlocations=String[],
    verbose::Bool=false)

## Input Parameters

- `c1`: Single Parameter

## Description

Cubic alpha `(α(T))` model. Default for `sCPA` EoS.
```
αᵢ = (1+c¹ᵢ(1-√(Trᵢ)))^2
```

"""
sCPAAlpha
default_locations(::Type{sCPAAlpha}) = ["SAFT/CPA/sCPA/sCPA_like.csv"]
