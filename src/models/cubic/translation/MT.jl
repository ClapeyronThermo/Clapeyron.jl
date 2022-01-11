abstract type MTTranslationModel <: TranslationModel end

struct MTTranslationParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple MTTranslation MTTranslationModel MTTranslationParam

export MTTranslation
function MTTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = MTTranslationParam(acentricfactor)
    model = MTTranslation(packagedparams, verbose=verbose)
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::MTTranslation)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    ω  = translation_model.params.acentricfactor.values
    c = zeros(typeof(T),length(Tc))
    for i ∈ @comps
        ωi = ω[i]
        Zc = evalpoly(ωi,(0.289,-0.0701,-0.0207))
        β  = -10.2447-28.6312*ωi
        Tci = Tc[i]
        RTp = R̄*Tci/Pc[i]
        t0 = RTp*evalpoly(ωi,(-0.014471,0.067498,-0.084852,0.067298,-0.017366))
        tc = RTp*(0.3074-Zc)
        Tr = T/Tci
        c[i] = t0+(tc-t0)*exp(β*abs(1-Tr))
    end
    return c
end