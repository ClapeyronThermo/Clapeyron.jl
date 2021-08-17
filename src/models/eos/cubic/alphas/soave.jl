abstract type SoaveAlphaModel <: AlphaModel end

struct SoaveAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

abstract type PRModel <: ABCubicModel end
@newmodelsimple SoaveAlpha SoaveAlphaModel SoaveAlphaParam

export SoaveAlpha
function SoaveAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = SoaveAlphaParam(acentricfactor)
    model = SoaveAlpha(packagedparams, verbose=verbose)
    return model
end

 