abstract type PenelouxTranslationModel <: TranslationModel end

struct PenelouxTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
    v_shift::SingleParam{Float64}
end

@newmodelsimple PenelouxTranslation PenelouxTranslationModel PenelouxTranslationParam

"""

    PenelouxTranslation <: PenelouxTranslationModel

    PenelouxTranslation(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`

## Model Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`


## Description

Peneloux Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
cᵢ = -0.252*RTcᵢ/Pcᵢ*(1.5448Zcᵢ - 0.4024)
Zcᵢ = Pcᵢ*Vcᵢ/(RTcᵢ)
```

## Model Construction Examples
```
# Using the default database
translation = PenelouxTranslation("water") #single input
translation = PenelouxTranslation(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
translation = PenelouxTranslation(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/Vc.csv"])

# Passing parameters directly
translation = PenelouxTranslation(["neon","hydrogen"];userlocations = (;Vc = [4.25e-5, 6.43e-5]))
```

## References

1. Péneloux A, Rauzy E, Fréze R. (1982) A consistent correction for Redlich‐Kwong‐Soave volumes. Fluid Phase Equilibria 1, 8(1), 7–23. [doi:10.1016/0378-3812(82)80002-2](https://doi.org/10.1016/0378-3812(82)80002-2)

"""
PenelouxTranslation

export PenelouxTranslation
default_locations(::Type{PenelouxTranslation}) = critical_data()
default_references(::Type{PenelouxTranslation}) = ["10.1016/0378-3812(82)80002-2"]
function transform_params(::Type{PenelouxTranslation},params,components)
    v_shift = SingleParam("Volume shift",components,zeros(length(components)))
    v_shift.ismissingvalues .= true
    params["v_shift"] = v_shift
    return params
end

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