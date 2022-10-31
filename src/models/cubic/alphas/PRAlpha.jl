PRAlpha_SETUP = ModelOptions(
        :PRAlpha;
        supertype=SoaveAlphaModel,
        locations=["properties/critical.csv"],
        inputparams=[
            ParamField(:w, SingleParam{Float64}),
        ],
        params=[
            ParamField(:acentricfactor, SingleParam{Float64}),
        ],
        mappings=[
            ModelMapping([:w], [:acentricfactor]),
        ],
    )

createmodel(PRAlpha_SETUP; verbose=true)
export PRAlpha

"""
    PRAlpha <: SoaveAlphaModel
    
    PRAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `w`: Single Parameter (`Float64`)

## Model Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`PR`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.37464 + 1.54226ωᵢ - 0.26992ωᵢ^2
```

"""
PRAlpha

@inline α_m(model,::PRAlpha) = (0.37464,1.54226,-0.26992)
