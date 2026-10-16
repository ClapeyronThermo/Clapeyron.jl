struct MC3PRAlphaParam{T} <: ParametricEoSParam{T}
    acentricfactor::SingleParam{T}
    c1::SingleParam{T}
    c2::SingleParam{T}
    c3::SingleParam{T}
end

MC3PRAlphaParam(ω,c1,c2,c3) = build_parametric_param(MC3PRAlphaParam,ω,c1,c2,c3)

@newmodelsimple MC3PRAlpha MathiasCopemanAlphaModel MC3PRAlphaParam

"""
    MC3PRAlpha <: MathiasCopemanAlphaModel

    MC3PRAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Predictive Mathias-Copeman alpha `(α(T))` model for PR models:
```
αᵢ = [1 + c₁ᵢ(1-√(Trᵢ)) + c₃ᵢ(1-√(Trᵢ))² + c₃ᵢ(1-√(Trᵢ))³)²
Trᵢ = T/Tcᵢ
c₁ᵢ = 0.396 + 1.3644*ωᵢ
c₃ᵢ = −0.0964 + 0.6593ωᵢ - 1.0793ωᵢ²
c₃ᵢ = 0.1656 − 0.0609ωᵢ + 1.1139ωᵢ²
```

For light gases (methane, nitrogen), custom parameters are used.

## Model Construction Examples
```
# Using the default database
alpha = MC3PRAlpha("water") #single input
alpha = MC3PRAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = MC3PRAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = MC3PRAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

## References

1.  P. M. Mathias and T. W. Copeman: Extension of the Peng-Robinson Equation of State to Complex Mixtures: Evaluation of the Various Forms of the Local Composition Concept, Fluid Phase Equilibria 13 (1983) 91–108, [doi:10.1016/0378-3812(83)80084-3](https://doi.org/10.1016/0378-3812(83)80084-3)
2.  Novak, N., Louli, V., & Voutsas, E. (2019). Prediction of Vapor–Liquid Equilibrium and Thermodynamic Properties of Natural Gas and Gas Condensates. Industrial & Engineering Chemistry Research, 58(17), 7370–7388. [doi:10.1021/acs.iecr.9b00756](https://doi.org/10.1021/acs.iecr.9b00756)
"""
MC3PRAlpha

default_locations(::Type{MC3PRAlpha}) = vcat(critical_data(),"alpha/MC3PRAlpha/MC3PRAlpha_like.csv")
default_ignore_missing_singleparams(::Type{MC3PRAlpha}) = ["Vc","c1","c2","c3","acentricfactor"]
default_references(::Type{MC3PRAlpha}) = ["10.1016/0378-3812(83)80084-3","10.1021/acs.iecr.9b00756"]
function transform_params(::Type{MC3PRAlpha},params,components)
    nc = length(components)
    ω = get!(params,"acentricfactor") do
        SingleParam("acentricfactor",components,zeros(nc),fill(true,nc))
    end

    c1 = get!(params,"c1") do
        SingleParam("c1",components,zeros(nc),fill(true,nc))
    end

    c2 = get!(params,"c2") do
        SingleParam("c2",components,zeros(nc),fill(true,nc))
    end

    c3 = get!(params,"c3") do
        SingleParam("c3",components,zeros(nc),fill(true,nc))
    end

    for i in 1:nc
        if c1.ismissingvalues[i]
            ω.ismissingvalues[i] && throw(error("MC3PRAlpha: cannot estimate c1: missing acentricfactor parameter"))
            c1[i] = 0.396 + 1.3644*ω[i]
        end

        if c2.ismissingvalues[i]
            ω.ismissingvalues[i] && throw(error("MC3PRAlpha: cannot estimate c2: missing acentricfactor parameter"))
            c2[i] = evalpoly(ω[i],(−0.0964,0.6593,-1.0793))
        end

        if c3.ismissingvalues[i]
            ω.ismissingvalues[i] && throw(error("MC3PRAlpha: cannot estimate c3: missing acentricfactor parameter"))
            c3[i] = evalpoly(ω[i],(0.1656,−0.0609,1.1139))
        end
    end
    return params
end

transform_params(::Type{MC3PRAlphaParam},params,components) = transform_params(MC3PRAlpha,params,components)

function recombine_impl!(model::MC3PRAlpha)
    ω,c1,c2,c3 = model.params.acentricfactor,model.params.c1,model.params.c2,model.params.c3
    for i in 1:length(model)
        if c1.ismissingvalues[i]
            ω.ismissingvalues[i] && throw(error("MC3PRAlpha: cannot estimate c1: missing acentricfactor parameter"))
            c1[i] = 0.396 + 1.3644*ω[i]
        end

        if c2.ismissingvalues[i]
            ω.ismissingvalues[i] && throw(error("MC3PRAlpha: cannot estimate c2: missing acentricfactor parameter"))
            c2[i] = evalpoly(ω[i],(−0.0964,0.6593,-1.0793))
        end

        if c3.ismissingvalues[i]
            ω.ismissingvalues[i] && throw(error("MC3PRAlpha: cannot estimate c3: missing acentricfactor parameter"))
            c3[i] = evalpoly(ω[i],(0.1656,−0.0609,1.1139))
        end
    end
    return model
end

export MC3PRAlpha