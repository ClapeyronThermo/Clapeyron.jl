abstract type PTVAlphaModel <: GeneralizedSuaveAlphaModel end

const PTVAlphaParam = SimpleAlphaParam

@newmodelsimple PTVAlpha PTVAlphaModel PTVAlphaParam
export PTVAlpha

"""
    PTVAlpha <: PTVAlphaModel
    
    PTVAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters


- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`PTV`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.46283 + 3.58230Zcᵢ*ωᵢ - 8.19417(Zcᵢ*ωᵢ)^2
```

## Model Construction Examples
```
# Using the default database
alpha = PTVAlpha("water") #single input
alpha = PTVAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = PTVAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = PTVAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

"""
PTVAlpha
default_locations(::Type{PTVAlpha}) = critical_data()

@inline function α_m(model::PTVModel,alpha_model::PTVAlphaModel,i)
    Tc = model.params.Tc.values[i]
    Pc = model.params.Pc.values[i]
    Vc = model.params.Vc.values[i]
    Zc = Vc*Pc/(R̄*Tc)
    coeff = (0.46283,3.58230,8.19417)
    ω = alpha_model.params.acentricfactor.values[i]
    return evalpoly(ω*Zc,coeff)
end