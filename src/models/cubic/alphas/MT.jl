struct MTAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple MTAlpha SoaveAlphaModel MTAlphaParam
export MTAlpha

"""
    MTAlpha <: MTAlphaModel
    
    MTAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
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

## References

1. Magoulas, K., & Tassios, D. (1990). Thermophysical properties of n-Alkanes from C1 to C20 and their prediction for higher ones. Fluid Phase Equilibria, 56, 119–140. [doi:10.1016/0378-3812(90)85098-u](https://doi.org/10.1016/0378-3812(90)85098-u)

"""
MTAlpha

function MTAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose,ignore_headers = ONLY_ACENTRICFACTOR)
    acentricfactor = params["acentricfactor"]
    packagedparams = MTAlphaParam(acentricfactor)
    model = MTAlpha(packagedparams, verbose=verbose)
    return model
end

@inline α_m(model,::MTAlpha) = (0.384401,1.52276,-0.213808,0.034616,-0.001976)
