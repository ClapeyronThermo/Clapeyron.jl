abstract type ConstantTranslationModel <: TranslationModel end

struct ConstantTranslationParam <: EoSParam
    v_shift::SingleParam{Float64}
end

@newmodelsimple ConstantTranslation ConstantTranslationModel ConstantTranslationParam

"""
    ConstantTranslation <: ConstantTranslationModel
    ConstantTranslation(components;
    userlocations = String[],
    verbose::Bool=false)
## Input Parameters

- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³·mol⁻¹]`

## Description
Constant Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
```
where `cᵢ` is constant. 
It does not have parameters by default, the volume shifts must be user-supplied.

## Model Construction Examples
```
# Using user-provided parameters

# Passing files or folders
translation = ConstantTranslation(["neon","hydrogen"]; userlocations = ["path/to/my/db","properties/critical.csv"])

# Passing parameters directly
translation = ConstantTranslation(["neon","hydrogen"];userlocations = (;Vc = [4.25e-5, 6.43e-5]))
```
"""
ConstantTranslation

export ConstantTranslation

recombine_translation!(model::CubicModel,translation_model::ConstantTranslation) = translation_model

function translation(model::CubicModel,V,T,z,translation_model::ConstantTranslationModel)
    return translation_model.params.v_shift.values
end