abstract type sCPAAlphaModel <: CPAAlphaModel end

@newmodelsimple sCPAAlpha sCPAAlphaModel CPAAlphaParam

export sCPAAlpha
function sCPAAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["SAFT/sCPA/sCPA_like.csv"]; userlocations=userlocations, ignore_missing_singleparams=["Mw"], verbose=verbose)
    c1 = params["c1"]
    packagedparams = CPAAlphaParam(c1)
    model = CPAAlpha(packagedparams, verbose=verbose)
    return model
end