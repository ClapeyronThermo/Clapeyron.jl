abstract type PenelouxTranslationModel <: TranslationModel end

struct PenelouxTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
end

@newmodelsimple PenelouxTranslation PenelouxTranslationModel PenelouxTranslationParam

export PenelouxTranslation
function PenelouxTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Vc = params["vc"]
    packagedparams = PenelouxTranslationParam(Vc)
    model = PenelouxTranslation(packagedparams, verbose=verbose)
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::PenelouxTranslation)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    c = zeros(eltype(Tc),length(Tc))
    for i ∈ @comps
        Tci = Tc[i]
        Pci = Pc[i]
        RT = Tci*R̄
        Zc = Pci*Vc[i]/RT
        c[i] = -0.252*RT/Pci*(1.5448*Zc-0.4024)
    end
    return c
end