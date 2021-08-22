abstract type RKAlphaModel <: AlphaModel end

struct RKAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

has_sites(::Type{<:AlphaModel})=false

@newmodelsimple RKAlpha RKAlphaModel RKAlphaParam

export RKAlpha
function RKAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = RKAlphaParam(acentricfactor)
    model = RKAlpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::CubicModel,V,T,z,alpha_model::RKAlphaModel)
    Tc = model.params.Tc.values
    α = @. 1 /√(T/Tc)
    return α
end