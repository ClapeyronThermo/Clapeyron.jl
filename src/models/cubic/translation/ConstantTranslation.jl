abstract type ConstantTranslationModel <: TranslationModel end

struct ConstantTranslationParam <: EoSParam
    v_shift::SingleParam{Float64}
end

@newmodelsimple ConstantTranslation ConstantTranslationModel ConstantTranslationParam

"""
    ConstantTranslation <: ConstantTranslationModel
    ConstantTranslation(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)
## Input Parameters

- `v_shift`: Single Parameter (`Float64`) - Volume shift `[m³/mol]`

## Description
Constant Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
```
where `cᵢ` is constant. 
It does not have parameters by default, the volume shifts must be user-supplied.
"""
ConstantTranslation

export ConstantTranslation

function ConstantTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, String[]; userlocations=userlocations, verbose=verbose)
    c = params["v_shift"]
    packagedparams = ConstantTranslationParam(c)
    model = ConstantTranslation(packagedparams, verbose=verbose)
    return model
end

recombine_translation!(model::CubicModel,translation_model::ConstantTranslation) = translation_model

function translation(model::CubicModel,V,T,z,translation_model::ConstantTranslationModel)
    return translation_model.params.v_shift.values
end
