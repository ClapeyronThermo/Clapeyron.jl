abstract type BMAlphaModel <: AlphaModel end

struct BMAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple BMAlpha BMAlphaModel BMAlphaParam
export BMAlpha

"""
    BMAlpha <: BMAlphaModel

    MTAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters


- `acentricfactor`: Single Parameter (`Float64`)

## Description

Boston Mathias Cubic alpha `(α(T))` model.
```
if Trᵢ > 1
    αᵢ = (exp((1-2/(2+mᵢ))*(1-Trᵢ^(1+mᵢ/2))))^2
else
    αᵢ = (1+mᵢ*(1-√(Trᵢ)))^2

Trᵢ = T/Tcᵢ

for PR models:
    mᵢ = 0.37464 + 1.54226ωᵢ - 0.26992ωᵢ^2
for RK models:
    mᵢ = 0.480 + 1.547ωᵢ - 0.176ωᵢ^2
```

## References

1. .M. Boston, P.M. Mathias, Proceedings of the 2nd International Conference on Phase Equilibria and Fluid Properties in the Chemical Process Industries, West Berlin, March, 1980, pp. 823–849

"""
BMAlpha

function BMAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose,ignore_headers = ONLY_ACENTRICFACTOR)
    acentricfactor = params["acentricfactor"]
    packagedparams = BMAlphaParam(acentricfactor)
    model = BMAlpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::RKModel,V,T,z,alpha_model::BMAlphaModel)
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        ωi = ω[i]
        m = evalpoly(ωi,(0.480,1.547,-0.176))
        Tr = T/Tc[i]
        α[i] = ifelse(Tr>1,(exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2,(1+m*(1-√(Tr)))^2)
    end
    return α
end

function α_function(model::PRModel,V,T,z,alpha_model::BMAlphaModel)
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    #m  = @. 0.37464+1.54226*ω-0.26992*ω^2
    #α  = @. (Tr>1)*(exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2+(Tr<=1)*(1+m*(1-√(Tr)))^2
    α = zeros(typeof(T),length(Tc))
    for i in @comps
        ωi = ω[i]
        m = evalpoly(ωi,(0.37464,1.54226,-0.26992))
        Tr = T/Tc[i]
        α[i] = ifelse(Tr>1,(exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2,(1+m*(1-√(Tr)))^2)
    end
    return α
end
