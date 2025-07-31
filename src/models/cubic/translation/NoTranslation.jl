abstract type NoTranslationModel <: TranslationModel end
@newmodelsingleton NoTranslation NoTranslationModel
export NoTranslation

"""
    NoTranslation(args...) <: TranslationModel

## Input Parameters

None

## Description

Default volume translation model for cubic models. It performs no translation:

```
V = V₀ + mixing_rule(cᵢ)
cᵢ = 0 ∀ i
```

## Model Construction Examples
```
# Because this model does not have parameters, all those constructors are equivalent:
translation = NoTranslation()
translation = NoTranslation("water")
translation = NoTranslation(["water","carbon dioxide"])
```
"""
NoTranslation

function translation(model::CubicModel,V,T,z,translation_model::NoTranslation)
    return FillArrays.Zeros{Float64}(length(z))
end

recombine_translation!(model::CubicModel,translation_model::NoTranslation) = translation_model