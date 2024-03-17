abstract type PTVAlphaModel <: AlphaModel end

const PTVAlphaParam = SimpleAlphaParam

@newmodelsimple PTVAlpha PTVAlphaModel PTVAlphaParam
export PTVAlpha

"""
    PTVAlpha <: PTVAlphaModel
    
    PTVAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters


- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`PTV`](@ref) EoS.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.46283 + 3.58230Zcᵢ*ωᵢ - 8.19417(Zcᵢ*ωᵢ)^2
```

## Model Construction Examples
```
# Using the default database
alpha = PTVAlpha("water") #single input
alpha = PTVAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = PTVAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = PTVAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

"""
PTVAlpha
default_locations(::Type{PTVAlpha}) = critical_data()

function α_function(model::CubicModel,V,T,z,alpha_model::PTVAlphaModel)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = model.params.Vc.values
    Zc = @. Vc*Pc/(R̄*Tc)
    ω  = alpha_model.params.acentricfactor.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        Zci = Vc[i]*Pc[i]/(R̄*Tc[i])
        coeff = (0.46283,3.58230*Zci,8.19417*Zci^2)
        ωi = ω[i]
        Tr = T/Tc[i]
        m = evalpoly(ωi,coeff)
        α[i] = (1+m*(1-√(Tr)))^2
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::PTVAlphaModel)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    Vc = model.params.Vc.values[1]
    Zc = Vc*Pc/(R̄*Tc)
    ω  = alpha_model.params.acentricfactor.values[1]
    coeff = (0.46283,3.58230*Zc,8.19417*Zc^2)
    Tr = T/Tc
    m = evalpoly(ω,coeff)
    α  = (1+m*(1-√(Tr)))^2
    return α
end