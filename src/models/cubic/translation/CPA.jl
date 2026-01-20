abstract type CPATranslationModel <: TranslationModel end

struct CPATranslationParam <: EoSParam
    c::SingleParam{Float64}
end

@newmodelsimple CPATranslation CPATranslationModel CPATranslationParam

export CPATranslation
function CPATranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["SAFT/sCPA/sCPA_like.csv"]; userlocations=userlocations, verbose=verbose)
    c = params["c"]
    c.values .*= 1e-6
    packagedparams = CPATranslationParam(c)
    model = CPATranslation(packagedparams, verbose=verbose)
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::CPATranslation)
    c = translation_model.params.c.values
    return c
end

function translation2(model::CubicModel,V,T,z,translation_model::CPATranslation,a,b,α)
    return dot(translation_model.params.c.values.values,z)
end