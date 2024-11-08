abstract type sCPAAlphaModel <: CPAAlphaModel end

@newmodelsimple sCPAAlpha sCPAAlphaModel CPAAlphaParam

"""
    sCPAAlpha <: sCPAAlphaModel
    
    sCPAAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `c1`: Single Parameter

## Description

Cubic alpha `(α(T))` model. Default for `sCPA` EoS.
```
αᵢ = (1+c¹ᵢ(1-√(Trᵢ)))^2
```

## Model Construction Examples
```
# Using the default database
alpha = sCPAAlpha("water") #single input
alpha = sCPAAlpha(["water","carbon dioxide"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = sCPAAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","scpa/alpha.csv"])

# Passing parameters directly
alpha = sCPAAlpha(["water","carbon dioxide"];userlocations = (;c1 = [0.67,0.76]))
```

"""
sCPAAlpha
default_locations(::Type{sCPAAlpha}) = ["SAFT/CPA/sCPA/sCPA_like.csv"]