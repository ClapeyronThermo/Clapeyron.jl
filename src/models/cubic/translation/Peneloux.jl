abstract type PenelouxTranslationModel <: TranslationModel end

struct PenelouxTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
    v_shift::SingleParam{Float64}
end

@newmodelsimple PenelouxTranslation PenelouxTranslationModel PenelouxTranslationParam

"""

    PenelouxTranslation <: PenelouxTranslationModel

    PenelouxTranslation(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³/mol]`

## Model Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³/mol]`
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³/mol]`


## Description

Peneloux Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
cᵢ = -0.252*RTcᵢ/Pcᵢ*(1.5448Zcᵢ - 0.4024)
Zcᵢ = Pcᵢ*Vcᵢ/(RTcᵢ)
```

## References

1. Péneloux A, Rauzy E, Fréze R. (1982) A consistent correction for Redlich‐Kwong‐Soave volumes. Fluid Phase Equilibria 1, 8(1), 7–23. [doi:10.1016/0378-3812(82)80002-2](https://doi.org/10.1016/0378-3812(82)80002-2)

"""
PenelouxTranslation

export PenelouxTranslation
function PenelouxTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose,ignore_headers = ONLY_VC)
    Vc = params["Vc"]
    c = SingleParam("Volume shift",components,zeros(length(components)))
    c.ismissingvalues .= true
    packagedparams = PenelouxTranslationParam(Vc,c)
    model = PenelouxTranslation(packagedparams, verbose=verbose)
    return model
end

doi(::PenelouxTranslation) = ["10.1016/0378-3812(82)80002-2"]

function translation(model::CubicModel,V,T,z,translation_model::PenelouxTranslation)
    c = translation_model.params.v_shift
    cmissing = c.ismissingvalues
    if any(cmissing)
        res = copy(c.values)
        translation!(model,V,T,z,translation_model,res)
    else
        res = c
    end
    return res
end

function recombine_translation!(model::CubicModel,translation_model::PenelouxTranslation)
    c = translation_model.params.v_shift
    translation!(model,0.0,0.0,0.0,translation_model,c.values)
    c.ismissingvalues .= false
    return translation_model
end

function translation!(model::CubicModel,V,T,z,translation_model::PenelouxTranslation,c)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    for i ∈ @comps
        Tci = Tc[i]
        Pci = Pc[i]
        RT = Tci*R̄
        Zc = Pci*Vc[i]/RT
        c[i] = -0.252*RT/Pci*(1.5448*Zc-0.4024)
    end
    return c
end
