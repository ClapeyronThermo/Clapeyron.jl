abstract type RackettTranslationModel <: TranslationModel end

struct RackettTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
    v_shift::SingleParam{Float64}
end

@newmodelsimple RackettTranslation RackettTranslationModel RackettTranslationParam

"""
    RackettTranslation <: RackettTranslationModel

    RackettTranslation(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`

## Model Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`

## Description

Rackett Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
cᵢ = 0.40768*RTcᵢ/Pcᵢ*(0.29441-Zcᵢ)
Zcᵢ = Pcᵢ*Vcᵢ/(RTcᵢ)
```

## Model Construction Examples
```
# Using the default database
translation = RackettTranslation("water") #single input
translation = RackettTranslation(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
translation = RackettTranslation(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/Vc.csv"])

# Passing parameters directly
translation = RackettTranslation(["neon","hydrogen"];userlocations = (;Vc = [4.25e-5, 6.43e-5]))
```

## References

1. Rackett, H. G. (1970). Equation of state for saturated liquids. Journal of Chemical and Engineering Data, 15(4), 514–517. [doi:10.1021/je60047a012](https://doi.org/10.1021/je60047a012)

"""
RackettTranslation

export RackettTranslation
default_locations(::Type{RackettTranslation}) = critical_data()
function transform_params(::Type{RackettTranslation},params,components)
    v_shift = SingleParam("Volume shift",components,zeros(length(components)))
    v_shift.ismissingvalues .= true
    params["v_shift"] = v_shift
    return params
end

doi(::RackettTranslation) = ["10.1016/0378-3812(82)80002-2"]

function translation!(model::CubicModel,V,T,z,translation_model::RackettTranslation,c)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    for i ∈ @comps
        Tci = Tc[i]
        Pci = Pc[i]
        RT = Tci*R̄
        Zc = Pci*Vc[i]/RT
        c[i] = 0.40768*RT/Pci*(0.29441-Zc)
    end
    return c
end

function translation(model::CubicModel,V,T,z,translation_model::RackettTranslation)
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

function recombine_translation!(model::CubicModel,translation_model::RackettTranslation)
    c = translation_model.params.v_shift
    translation!(model,0.0,0.0,0.0,translation_model,c.values)
    c.ismissingvalues .= false
    return translation_model
end