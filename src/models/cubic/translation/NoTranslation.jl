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

function NoTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    model = NoTranslation(NoTranslationParam())
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::NoTranslation)
    return FillArrays.Zeros{Float64}(length(z))
end

is_splittable(::NoTranslation) = false