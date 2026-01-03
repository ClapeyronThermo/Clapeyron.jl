abstract type VTPRTranslationModel <: RackettTranslationModel end

@newmodelsimple VTPRTranslation RackettTranslationModel RackettTranslationParam

"""
    VTPRTranslation <: RackettTranslationModel

    VTPRTranslation(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`

## Model Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`
- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`

## Description

VTPR Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
cᵢ = 0.252*RTcᵢ/Pcᵢ*(1.5448Zcᵢ - 0.4024)
Zcᵢ = Pcᵢ*Vcᵢ/(RTcᵢ)
```

## Model Construction Examples
```
# Using the default database
translation = VTPRTranslation("water") #single input
translation = VTPRTranslation(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
translation = VTPRTranslation(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/Vc.csv"])

# Passing parameters directly
translation = VTPRTranslation(["neon","hydrogen"];userlocations = (;Vc = [4.25e-5, 6.43e-5]))
```

## References

1. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4) 

"""
VTPRTranslation

export VTPRTranslation
default_locations(::Type{VTPRTranslation}) = vcat("cubic/VTPR/VTPRTranslation_like.csv",critical_data())
default_ignore_missing_singleparams(::Type{VTPRTranslation}) = ["v_shift,Vc"]

function transform_params(::Type{VTPRTranslation},params,components)
    v_shift = get!(params,"v_shift") do
        SingleParam("Volume shift",components,zeros(length(components)))
    end
    params["v_shift"] = v_shift
    return params
end

doi(::VTPRTranslation) = ["10.1016/s0378-3812(01)00626-4"]

function translation!(model::CubicModel,V,T,z,translation_model::VTPRTranslation,c)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    missing_c = translation_model.params.v_shift.ismissingvalues
    missing_vc = translation_model.params.Vc.ismissingvalues
    for i ∈ @comps
        if missing_c[i] && !missing_vc[i]
            Tci = Tc[i]
            Pci = Pc[i]
            RT = Tci*R̄
            Zc = Pci*Vc[i]/RT
            c[i] = 0.252*RT/Pci*(1.5448Zc - 0.4024)
        elseif missing_vc[i] && missing_c[i]
            throw(MissingException("VTPRTranslation: cannot estimate v_shift: missing Vc parameter."))
        end
    end
    return c
end

function translation(model::CubicModel,V,T,z,translation_model::VTPRTranslation)
    c = translation_model.params.v_shift
    return c
end

function recombine_translation!(model::CubicModel,translation_model::VTPRTranslation)
    c = translation_model.params.v_shift
    translation!(model,0.0,0.0,0.0,translation_model,c.values)
    c.ismissingvalues .= false
    return translation_model
end