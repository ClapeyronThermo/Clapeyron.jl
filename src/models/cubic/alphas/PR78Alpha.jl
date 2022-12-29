abstract type PR78AlphaModel <: AlphaModel end

struct PR78AlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple PR78Alpha PR78AlphaModel PR78AlphaParam
export PR78Alpha

"""
    PR78Alpha <: PR78AlphaModel
    
    PR78Alpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters
- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`PR78`](@ref) and [`EPPR78`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
if ωᵢ ≤ 0.491
    mᵢ = 0.37464 + 1.54226ωᵢ - 0.26992ωᵢ^2
else
    mᵢ = 0.379642 + 1.487503ωᵢ - 0.164423ωᵢ^2 - 0.016666ωᵢ^3
```
"""
PR78Alpha

function PR78Alpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose,ignore_headers = ONLY_ACENTRICFACTOR)
    acentricfactor = params["acentricfactor"]
    packagedparams = PR78AlphaParam(acentricfactor)
    model = PR78Alpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::CubicModel,V,T,z,alpha_model::PR78AlphaModel)
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        ωi = ω[i]
        m = ifelse(ωi<=0.491,
            evalpoly(ωi,(0.37464,1.54226,-0.26992)),
            evalpoly(ωi,(0.379642,1.487503,-0.164423,-0.016666)))
        Tr = T/Tc[i]
        α[i] = (1+m*(1-√(Tr)))^2
    end
    
    return α

end