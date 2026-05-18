abstract type BaledTranslationModel <: TranslationModel end

struct BaledTranslationParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple BaledTranslation BaledTranslationModel SimpleAlphaParam
export BaledTranslation

"""

BaledTranslation <: BaledTranslationModel

    BaledTranslationModel(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Baled Translation model for cubics:

```
V = V₀ + mixing_rule(cᵢ)
cᵢ = A + BTrᵢ
Aᵢ = k₀ + k₁*exp(-1/k₂fᵢ) + k₃*exp(-1/k₄fᵢ) + k₅*exp(-1/k₆fᵢ)
Bᵢ = j₀ + j₁*exp(-1/j₂fᵢ) + j₃*exp(-1/j₄fᵢ) + j₅*exp(-1/j₆fᵢ)
fᵢ = Mwᵢ*ωᵢ
```

## Model Construction Examples
```
# Using the default database
translation = BaledTranslation("water") #single input
translation = BaledTranslation(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
translation = BaledTranslation(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
translation = BaledTranslation(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```


## References

1. Baled, H., Enick, R. M., Wu, Y., McHugh, M. A., Burgess, W., Tapriyal, D., & Morreale, B. D. (2012). Prediction of hydrocarbon densities at extreme conditions using volume-translated SRK and PR equations of state fit to high temperature, high pressure PVT data. Fluid Phase Equilibria, 317, 65–76. [doi:10.1016/j.fluid.2011.12.027](https://doi.org/10.1016/j.fluid.2011.12.027)

"""
BaledTranslation
default_locations(::Type{BaledTranslation}) = critical_data()
default_references(::Type{BaledTranslation}) = ["10.1016/j.fluid.2011.12.027"]

function translation(model::CubicModel,V,T,z,translation_model::BaledTranslation)
    Tc = model.params.Tc.values
    ω  = translation_model.params.acentricfactor.values
    MW = model.params.Mw.values
    k0A,k0B,kA,kB = baled_translations_consts(model)
    c = zeros(Base.promote_eltype(translation_model,Tc,T),length(Tc))
    for i ∈ @comps
        Tci = Tc[i]
        Tri = T/Tci
        fi = -1/(ω[i]*1000*MW[i])
        A = k0A + kA[1]*exp(fi/kA[2]) + kA[3]*exp(fi/kA[4]) + kA[5]*exp(fi/kA[6])
        B = k0B + kB[1]*exp(fi/kB[2]) + kB[3]*exp(fi/kB[4]) + kB[5]*exp(fi/kB[6])
        c[i] = (A + B*Tri)*1e-6
    end
    return c
end

function translation2(model::CubicModel,V,T,z,translation_model::BaledTranslation,a,b,α)
    Tc = model.params.Tc.values
    ω  = translation_model.params.acentricfactor.values
    MW = model.params.Mw.values
    k0A,k0B,kA,kB = baled_translations_consts(model)
    c = zero(Base.promote_eltype(model,T,z))
    for i ∈ @comps
        Tci = Tc[i]
        Tri = T/Tci
        fi = -1/(ω[i]*1000*MW[i])
        A = k0A + kA[1]*exp(fi/kA[2]) + kA[3]*exp(fi/kA[4]) + kA[5]*exp(fi/kA[6])
        B = k0B + kB[1]*exp(fi/kB[2]) + kB[3]*exp(fi/kB[4]) + kB[5]*exp(fi/kB[6])
        c += z[i]*(A + B*Tri)
    end
    return c*1e-6
end

function baled_translations_consts(model::PRModel)
    k0A = −4.1034
    kA = (31.723, 0.0531, 188.68, 0.0057, 20.196, 0.0003)
    k0B = −0.3489
    kB = (−28.547, 0.0687, −817.73,0.0007, −65.067, 0.0076)
    return k0A,k0B,kA,kB
end

function baled_translations_consts(model::RKModel)
    k0A = 0.23
    kA = (46.843, 0.0571, 23.161, 0.0003, 267.40, 0.0053)
    k0B = −0.3471
    kB = (−29.748, 0.0644, −347.04, 0.0010, −88.547, 0.0048)
    return k0A,k0B,kA,kB
end