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