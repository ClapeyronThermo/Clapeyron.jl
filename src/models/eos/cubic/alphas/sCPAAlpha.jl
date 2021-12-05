abstract type sCPAAlphaModel <: AlphaModel end

struct sCPAAlphaParam <: EoSParam
    c1::SingleParam{Float64}
end

@newmodelsimple sCPAAlpha sCPAAlphaModel sCPAAlphaParam

export sCPAAlpha
function sCPAAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["SAFT/CPA/sCPA/sCPA_like.csv"]; userlocations=userlocations, ignore_missing_singleparams=["Mw"], verbose=verbose)
    c1 = params["c1"]
    packagedparams = sCPAAlphaParam(c1)
    model = sCPAAlpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::CubicModel,V,T,z,alpha_model::sCPAAlphaModel)
    Tc = model.params.Tc.values
    Tr = @. T/Tc
    c1  = alpha_model.params.c1.values
    α  = @. (1+c1*(1-√(Tr)))^2
    return α
end