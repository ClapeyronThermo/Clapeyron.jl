const PatelTejaAlphaParam = SimpleAlphaParam

@newmodelsimple PatelTejaAlpha SoaveAlphaModel PatelTejaAlphaParam
export PatelTejaAlpha

"""
    PatelTejaAlpha <: SoaveAlphaModel
    
    PatelTejaAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`PatelTeja`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.452413 + 1.30982ωᵢ - 0.295937ωᵢ^2
```

## Model Construction Examples
```
# Using the default database
alpha = PatelTejaAlpha("water") #single input
alpha = PatelTejaAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = PatelTejaAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = PatelTejaAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

"""
PatelTejaAlpha
default_locations(::Type{PatelTejaAlpha}) = critical_data()

@inline α_m(model,::PatelTejaAlpha) = (0.452413,1.30982,-0.295937)
