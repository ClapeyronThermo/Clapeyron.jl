abstract type NoTranslationModel <: TranslationModel end

struct NoTranslation <: NoTranslationModel end

is_splittable(::NoTranslation) = false

export NoTranslation

"""
    NoTranslation(args...) <: TranslationModel

## Input Parameters

None

## Description

Default volume translation model for cubic models. it performs no translation:

```
V = V₀ + mixing_rule(cᵢ)
cᵢ = 0 ∀ i
```

"""
function NoTranslation(components;userlocations = String[],verbose = false)
    return NoTranslation()
end

function translation(model::CubicModel,V,T,z,translation_model::NoTranslation)
    return FillArrays.Zeros{Float64}(length(z))
end

recombine_translation!(model::CubicModel,translation_model::NoTranslation) = translation_model
