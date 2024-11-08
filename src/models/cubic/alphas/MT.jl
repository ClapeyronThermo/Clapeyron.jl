const MTAlphaParam = SimpleAlphaParam
@newmodelsimple MTAlpha SoaveAlphaModel MTAlphaParam
export MTAlpha

"""
    MTAlpha <: MTAlphaModel
    
    MTAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Magoulas & Tassios Cubic alpha `(α(T))` model. Default for [`UMRPR`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.384401 + 1.52276ωᵢ - 0.213808ωᵢ^2 + 0.034616ωᵢ^3 - 0.001976ωᵢ^4 
```

## Model Construction Examples
```
# Using the default database
alpha = MTAlpha("water") #single input
alpha = MTAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = MTAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = MTAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

## References

1. Magoulas, K., & Tassios, D. (1990). Thermophysical properties of n-Alkanes from C1 to C20 and their prediction for higher ones. Fluid Phase Equilibria, 56, 119–140. [doi:10.1016/0378-3812(90)85098-u](https://doi.org/10.1016/0378-3812(90)85098-u)

"""
MTAlpha
default_locations(::Type{MTAlpha}) = critical_data()

@inline α_m(model,::MTAlpha) = (0.384401,1.52276,-0.213808,0.034616,-0.001976)