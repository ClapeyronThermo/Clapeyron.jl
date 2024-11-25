abstract type CPAAlphaModel <: AlphaModel end

struct CPAAlphaParam <: EoSParam
    c1::SingleParam{Float64}
end

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

function α_function!(α,model::CubicModel,alpha_model::CPAAlphaModel,T)
    Tc = model.params.Tc.values
    c1  = alpha_model.params.c1.values
    for i in @comps
        Tr = T/Tc[i]
        α[i] = (1+c1[i]*(1-√(Tr)))^2
    end
    return α
end

function α_function(model::CubicModel,alpha_model::CPAAlphaModel,T,i::Int)
    Tc = model.params.Tc.values[i]
    c1  = alpha_model.params.c1.values[i]
    Tr = T/Tc
    α = (1+c1*(1-√(Tr)))^2
end