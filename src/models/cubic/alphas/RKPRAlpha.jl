abstract type RKPRAlphaModel <: AlphaModel end

struct RKPRAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple RKPRAlpha RKPRAlphaModel RKPRAlphaParam
export RKPRAlpha

"""
    RKPRAlpha <: RKPRAlphaModel
    
    RKPRAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters


- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`RKPR`](@ref) EoS.
```
αᵢ = (3/(2 + Trᵢ))^kᵢ
kᵢ = (12.504Zc -2.7238) + (7.4513*Zc + 1.9681)ωᵢ + (-2.4407*Zc + 0.0017)ωᵢ^2
Trᵢ = T/Tcᵢ
```
"""
RKPRAlpha

function RKPRAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose,ignore_headers = ONLY_ACENTRICFACTOR)
    acentricfactor = params["acentricfactor"]
    packagedparams = RKPRAlphaParam(acentricfactor)
    model = RKPRAlpha(packagedparams, verbose=verbose)
    return model
end
#ideally, k should be fitted to Tr = 0.7
function α_function(model::CubicModel,V,T,z,alpha_model::RKPRAlphaModel)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = model.params.Vc.values
    c = model.params.c.values
    ω  = alpha_model.params.acentricfactor.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        δ = c[i]
        d = (1 + δ*δ)/(1+δ)
        y = 1 + cbrt(2*(1+δ)) + cbrt(4/(1+δ))
        Zc = y/(3y + d - 1)
        ωi = ω[i]
        Tr = T/Tc[i]
        ki = evalpoly(ωi,(12.504*Zc -2.7238,7.4513*Zc + 1.9681,-2.4407*Zc + 0.0017))
        α[i] = (3/(2 + Tr))^ki
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::RKPRAlphaModel)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    ω  = alpha_model.params.acentricfactor.values[1]
    δ = model.params.c.values[1]
    d = (1 + δ*δ)/(1+δ)
    y = 1 + cbrt(2*(1+δ)) + cbrt(4/(1+δ))
    Zc = y/(3y + d - 1)
    Tr = T/Tc
    ki = evalpoly(ω,(12.504*Zc -2.7238,7.4513*Zc + 1.9681,-2.4407*Zc + 0.0017))
    return (3/(2 + Tr))^ki
end

