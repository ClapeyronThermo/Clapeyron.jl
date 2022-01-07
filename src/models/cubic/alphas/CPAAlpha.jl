abstract type CPAAlphaModel <: AlphaModel end

struct CPAAlphaParam <: EoSParam
    c1::SingleParam{Float64}
end

@newmodelsimple CPAAlpha CPAAlphaModel CPAAlphaParam

export CPAAlpha
function CPAAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["SAFT/CPA/CPA_like.csv"]; userlocations=userlocations, ignore_missing_singleparams=["Mw"], verbose=verbose)
    c1 = params["c1"]
    packagedparams = CPAAlphaParam(c1)
    model = CPAAlpha(packagedparams, verbose=verbose)
    return model
end

