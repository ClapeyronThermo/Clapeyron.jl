abstract type CPAAlphaModel <: GeneralizedSuaveAlphaModel end

#=
struct CPAAlphaParam <: EoSParam
    c1::SingleParam{Float64}
end
moved to alphas.jl
=#

@newmodelsimple CPAAlpha CPAAlphaModel CPAAlphaParam
export CPAAlpha

"""
    CPAAlpha <: CPAAlphaModel
    
    CPAAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `c1`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for `CPA` EoS.
```
αᵢ = (1+c¹ᵢ(1-√(Trᵢ)))^2
```

## Model Construction Examples
```
# Using the default database
alpha = CPAAlpha("water") #single input
alpha = CPAAlpha(["water","carbon dioxide"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = CPAAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","cpa/alpha.csv"])

# Passing parameters directly
alpha = CPAAlpha(["water","carbon dioxide"];userlocations = (;c1 = [0.67,0.76]))
```
"""
CPAAlpha
default_locations(::Type{CPAAlpha}) = ["SAFT/CPA/CPA_like.csv"]

@inline function α_m(model::CubicModel,alpha_model::CPAAlphaModel,i)
    return alpha_model.params.c1.values[i]
end

const sCPAAlphaModel = CPAAlphaModel

@newmodelsimple sCPAAlpha CPAAlphaModel CPAAlphaParam

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

@newmodelsimple eCPAAlpha CPAAlphaModel CPAAlphaParam

"""
    eCPAAlpha <: CPAAlphaModel
    
    eCPAAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `c1`: Single Parameter

## Description

Cubic alpha `(α(T))` model. Default for `eCPA` EoS.
```
αᵢ = (1+c¹ᵢ(1-√(Trᵢ)))^2
```

## Model Construction Examples
```
# Using the default database
alpha = eCPAAlpha("water") #single input
alpha = eCPAAlpha(["water","carbon dioxide"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = eCPAAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","scpa/alpha.csv"])

# Passing parameters directly
alpha = eCPAAlpha(["water","carbon dioxide"];userlocations = (;c1 = [0.67,0.76]))
```

"""
eCPAAlpha
default_locations(::Type{eCPAAlpha}) = ["SAFT/CPA/eCPA/eCPA_like.csv"]