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

    Zc = @. 0.289-0.0701*ω-0.0207*ω^2
    β  = @. -10.2447-28.6312*ω

    t0 = @. R̄*Tc/Pc*(-0.014471+0.067498*ω-0.084852*ω^2+0.067298*ω^3-0.017366*ω^4)
    tc = @. R̄*Tc/Pc*(0.3074-Zc)

    Tr = @. T/Tc
    return @. t0+(tc-t0)*exp(β*abs(1-Tr))
end