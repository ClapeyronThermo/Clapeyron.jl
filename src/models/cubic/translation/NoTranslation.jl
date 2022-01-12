abstract type NoTranslationModel <: TranslationModel end

struct NoTranslationParam <: EoSParam
end

@newmodelsimple NoTranslation NoTranslationModel NoTranslationParam

export NoTranslation
function NoTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    model = NoTranslation(NoTranslationParam())
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::NoTranslation)
    return FillArrays.Zeros{Float64}(length(z))
end

is_splittable(::NoTranslation) = false