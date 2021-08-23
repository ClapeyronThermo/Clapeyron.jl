abstract type PR78AlphaModel <: AlphaModel end

struct PR78AlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple PR78Alpha PR78AlphaModel PR78AlphaParam

export PR78Alpha
function PR78Alpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = PR78AlphaParam(acentricfactor)
    model = PR78Alpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::CubicModel,V,T,z,alpha_model::PR78AlphaModel)
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    m = @. (ω<=0.491)*(0.37464+1.54226*ω-0.26992*ω^2)+(ω>0.491)*(0.379642+1.487503*ω-0.164423*ω^2-0.016666*ω^3)
    α = @. (1+m*(1-√(T/Tc)))^2
    return α
end