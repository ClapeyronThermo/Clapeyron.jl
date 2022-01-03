abstract type TwuAlphaModel <: AlphaModel end

struct TwuAlphaParam <: EoSParam
    M::SingleParam{Float64}
    N::SingleParam{Float64}
    L::SingleParam{Float64}
end

@newmodelsimple TwuAlpha TwuAlphaModel TwuAlphaParam

export TwuAlpha
function TwuAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv","alpha/Twu/Twu_like.csv"]; userlocations=userlocations, verbose=verbose)
    M = params["M"]
    N = params["N"]
    L = params["L"]
    packagedparams = TwuAlphaParam(M,N,L)
    model = TwuAlpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::CubicModel,V,T,z,alpha_model::TwuAlphaModel)
    Tc = model.params.Tc.values
    M  = alpha_model.params.M.values
    N  = alpha_model.params.N.values
    L  = alpha_model.params.L.values
    Tr = T ./ Tc
    α = @. Tr^(N*(M-1))*exp(L*(1-Tr^(N*M)))
    return α
end