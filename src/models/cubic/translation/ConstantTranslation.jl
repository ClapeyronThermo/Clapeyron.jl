abstract type ConstantTranslationModel <: TranslationModel end

ConstantTranslation_SETUP = ModelOptions(
        :ConstantTranslation;
        supertype=ConstantTranslationModel,
        params=[ParamField(:c, SingleParam{Float64})],
    )

createmodel(ConstantTranslation_SETUP; verbose=true)
export ConstantTranslation

"""

    ConstantTranslation <: ConstantTranslationModel

    ConstantTranslation(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `c`: Single Parameter (`Float64`) - Volume shift `[m³/mol]`

## Description

Constant Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
```
where `cᵢ` is constant. 

It does not have parameters by default, the volume shifts must be user-supplied.

"""
ConstantTranslation

function translation(model::CubicModel,V,T,z,translation_model::ConstantTranslationModel)
    return translation_model.params.c.values
end

