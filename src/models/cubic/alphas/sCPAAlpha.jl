abstract type sCPAAlphaModel <: CPAAlphaModel end


@newmodelsimple sCPAAlpha sCPAAlphaModel CPAAlphaParam

"""
    sCPAAlpha <: sCPAAlphaModel
    
    CPAAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `c1`: Single Parameter

## Description

Cubic alpha `(α(T))` model. Default for `sCPA` EoS.
```
αᵢ = (1+c¹ᵢ(1-√(Trᵢ)))^2
```

"""
sCPAAlpha

export sCPAAlpha
function sCPAAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false, kwargs...)
    params = getparams(components, ["SAFT/CPA/sCPA/sCPA_like.csv"]; userlocations=userlocations, ignore_missing_singleparams=["Mw"], verbose=verbose)
    c1 = params["c1"]
    packagedparams = CPAAlphaParam(c1)
    model = sCPAAlpha(packagedparams, verbose=verbose)
    return model
end
