abstract type SoaveAlphaModel <: AlphaModel end

const SoaveAlphaParam = SimpleAlphaParam

@newmodelsimple SoaveAlpha SoaveAlphaModel SoaveAlphaParam
export SoaveAlpha

"""
    SoaveAlpha <: SoaveAlphaModel
    
    SoaveAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters


- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`SRK`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.480 + 1.547ωᵢ - 0.176ωᵢ^2
```
To use different polynomial coefficients for `mᵢ`, overload `Clapeyron.α_m(::CubicModel,::SoaveAlphaModel) = (c₁,c₂,...cₙ)`

## Model Construction Examples
```
# Using the default database
alpha = SoaveAlpha("water") #single input
alpha = SoaveAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = SoaveAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = SoaveAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

"""
SoaveAlpha
default_locations(::Type{SoaveAlpha}) = critical_data()

@inline α_m(model::RKModel,::SoaveAlpha) = (0.480,1.547,-0.176)
@inline α_m(model::PRModel,::SoaveAlpha) = (0.37464,1.54226,-0.26992) #equal to PRAlpha
@inline α_m(model::vdWModel,::SoaveAlpha) = (0.4998,1.5928,0.19563,0.025)

function α_function(model::CubicModel,V,T,z,alpha_model::SoaveAlphaModel)
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    α = zeros(typeof(T*1.0),length(Tc))
    coeff = α_m(model,alpha_model)
    for i in @comps
        ωi = ω[i]
        Tr = T/Tc[i]
        m = evalpoly(ωi,coeff)
        α[i] = (1+m*(1-√(Tr)))^2
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::SoaveAlphaModel)
    Tc = model.params.Tc.values[1]
    ω  = alpha_model.params.acentricfactor.values[1]
    coeff = α_m(model,alpha_model)
    Tr = T/Tc
    m = evalpoly(ω,coeff)
    α  = (1+m*(1-√(Tr)))^2
    return α
end

const SRKModel = RK{I,SoaveAlpha,M,T} where {I,M,T}