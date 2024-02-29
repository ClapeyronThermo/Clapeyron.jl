const PRAlphaParam = SimpleAlphaParam

@newmodelsimple PRAlpha SoaveAlphaModel PRAlphaParam
export PRAlpha

"""
    PRAlpha <: SoaveAlphaModel
    
    PRAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description
Cubic alpha `(α(T))` model. Default for [`PR`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.37464 + 1.54226ωᵢ - 0.26992ωᵢ^2
```

It is equivalent to `SoaveAlpha`.

## Model Construction Examples
```
# Using the default database
alpha = PRAlpha("water") #single input
alpha = PRAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = PRAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = PRAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

"""
PRAlpha
default_locations(::Type{PRAlpha}) = critical_data()

@inline α_m(model,::PRAlpha) = (0.37464,1.54226,-0.26992)