abstract type BMAlphaModel <: AlphaModel end

const BMAlphaParam = SimpleAlphaParam

@newmodelsimple BMAlpha BMAlphaModel SimpleAlphaParam
export BMAlpha

"""
    BMAlpha <: BMAlphaModel

    BMAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Boston Mathias Cubic alpha `(α(T))` model. Identical to a Soave model when Tr <= 1, while having the correct behaviour at supercritical temperatures (via extrapolation at the critical point).
```julia
if Trᵢ > 1
    αᵢ = (exp((1-2/(2+mᵢ))*(1-Trᵢ^(1+mᵢ/2))))^2
else
    αᵢ = (1+mᵢ*(1-√(Trᵢ)))^2

Trᵢ = T/Tcᵢ

for PR models, RK, vdW models:
    mᵢ = m(Soave)
other cubic models:
    mᵢ = m(Leivobici)
```

## Model Construction Examples
```
# Using the default database
alpha = BMAlpha("water") #single input
alpha = BMAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = BMAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = BMAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

## References

1. .M. Boston, P.M. Mathias, Proceedings of the 2nd International Conference on Phase Equilibria and Fluid Properties in the Chemical Process Industries, West Berlin, March, 1980, pp. 823–849

"""
BMAlpha
default_locations(::Type{BMAlpha}) = critical_data()
default_ignore_missing_singleparams(::Type{BMAlpha}) = ["Vc"]

function α_function(model::DeltaCubicModel,V,T,z,alpha_model::BMAlphaModel)
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        ωi = ω[i]
        poly = α_m_BM(model,alpha_model,i)
        m = evalpoly(ωi,poly)
        Tr = T/Tc[i]
        if Tr > 1
            α[i] = (exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2
        else
            α[i] = (1+m*(1-sqrt(Tr)))^2
        end
    end
    return α
end

function α_function(model::DeltaCubicModel,V,T,z::SingleComp,alpha_model::BMAlphaModel)
    Tc = model.params.Tc.values[1]
    ωi  = alpha_model.params.acentricfactor.values[1]
    poly = α_m_BM(model,alpha_model,1)
    m = evalpoly(ωi,poly)
    Tr = T/Tc
    if Tr > 1
        α = (exp((1-2/(2+m))*(1-Tr^(1+m/2))))^2
    else
        α = (1+m*(1-sqrt(Tr)))^2
    end
    return α
end

@inline α_m_BM(model::RKModel,::BMAlphaModel,i) = (0.480,1.547,-0.176)
@inline α_m_BM(model::PRModel,::BMAlphaModel,i) = (0.37464,1.54226,-0.26992)
@inline α_m_BM(model::vdWModel,::BMAlphaModel,i) = (0.4998,1.5928,0.19563,0.025)
α_m_BM(model::DeltaCubicModel,::BMAlphaModel,i) = α_m_leibovici(model,i)