abstract type PR78AlphaModel <: GeneralizedSuaveAlphaModel end

const PR78AlphaParam = SimpleAlphaParam

@newmodelsimple PR78Alpha PR78AlphaModel PR78AlphaParam
export PR78Alpha

"""
    PR78Alpha <: PR78AlphaModel
    
    PR78Alpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters
- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`PR78`](@ref) and [`EPPR78`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
if ωᵢ ≤ 0.491
    mᵢ = 0.37464 + 1.54226ωᵢ - 0.26992ωᵢ^2
else
    mᵢ = 0.379642 + 1.487503ωᵢ - 0.164423ωᵢ^2 - 0.016666ωᵢ^3
```

## Model Construction Examples
```
# Using the default database
alpha = PR78Alpha("water") #single input
alpha = PR78Alpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = PR78Alpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = PR78Alpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```
"""
PR78Alpha
default_locations(::Type{PR78Alpha}) = critical_data()

@inline function α_m(model::PRModel,alpha_model::PR78AlphaModel,i)
    ωi  = alpha_model.params.acentricfactor.values[i]
    if ωi <= 0.491
        return evalpoly(ωi,(0.37464,1.54226,-0.26992))
    else
        return evalpoly(ωi,(0.379642,1.487503,-0.164423,-0.016666))
    end
end
