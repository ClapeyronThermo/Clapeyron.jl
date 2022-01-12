abstract type RackettTranslationModel <: TranslationModel end

struct RackettTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
end

@newmodelsimple RackettTranslation RackettTranslationModel RackettTranslationParam

export RackettTranslation
function RackettTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Vc = params["vc"]
    packagedparams = RackettTranslationParam(Vc)
    model = RackettTranslation(packagedparams, verbose=verbose)
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::RackettTranslation)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    c = zeros(eltype(Tc),length(Tc))
    for i ∈ @comps
        Tci = Tc[i]
        Pci = Pc[i]
        RT = Tci*R̄
        Zc = Pci*Vc[i]/RT
        c[i] = 0.40768*RT/Pci*(0.29441-Zc)
    end
    return c
end