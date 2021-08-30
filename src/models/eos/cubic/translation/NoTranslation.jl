abstract type NoTranslationModel <: TranslationModel end

struct NoTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
end

@newmodelsimple NoTranslation NoTranslationModel NoTranslationParam

export NoTranslation
function NoTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Vc = params["vc"]
    packagedparams = NoTranslationParam(Vc)
    model = NoTranslation(packagedparams, verbose=verbose)
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::NoTranslation)
    return zero(z)
end