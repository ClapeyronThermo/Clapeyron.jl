const PatelTejaAlphaParam = SimpleAlphaParam

@newmodelsimple PatelTejaAlpha SoaveAlphaModel PatelTejaAlphaParam
export PatelTejaAlpha

"""
    PatelTejaAlpha <: SoaveAlphaModel
    
    PatelTejaAlpha(components;
    userlocations=String[],
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

"""
PatelTejaAlpha
default_locations(::Type{PatelTejaAlpha}) = critical_data()


@inline α_m(model,::PatelTejaAlpha) = (0.452413,1.30982,-0.295937)