abstract type NoTranslationModel <: TranslationModel end

struct NoTranslationParam <: EoSParam
end

@newmodelsimple NoTranslation NoTranslationModel NoTranslationParam

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
NoTranslation

function NoTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false, kwargs...)
    model = NoTranslation(NoTranslationParam())
    return model
end

NoTranslation() = NoTranslation(NoTranslationParam())

function translation(model::CubicModel,V,T,z,translation_model::NoTranslation)
    return FillArrays.Zeros{Float64}(length(z))
end

recombine_translation!(model::CubicModel,translation_model::NoTranslation) = translation_model


is_splittable(::NoTranslation) = false
