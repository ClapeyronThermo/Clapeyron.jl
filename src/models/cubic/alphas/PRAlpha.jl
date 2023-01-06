struct PRAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple PRAlpha SoaveAlphaModel PRAlphaParam
export PRAlpha

"""
    PRAlpha <: SoaveAlphaModel
    
    PRAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
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
"""
PRAlpha

function PRAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose,ignore_headers = ONLY_ACENTRICFACTOR)
    acentricfactor = params["acentricfactor"]
    packagedparams = PRAlphaParam(acentricfactor)
    model = PRAlpha(packagedparams, verbose=verbose)
    return model
end

@inline α_m(model,::PRAlpha) = (0.37464,1.54226,-0.26992)