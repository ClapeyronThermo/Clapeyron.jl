struct PatelTejaAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple PatelTejaAlpha SoaveAlphaModel PatelTejaAlphaParam
export PatelTejaAlpha

"""
    PatelTejaAlpha <: SoaveAlphaModel
    
    PatelTejaAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `w`: Single Parameter (`Float64`)

## Model Parameters

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

function PatelTejaAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = PatelTejaAlphaParam(acentricfactor)
    model = PatelTejaAlpha(packagedparams, verbose=verbose)
    return model
end

@inline α_m(model,::PatelTejaAlpha) = (0.452413,1.30982,-0.295937)