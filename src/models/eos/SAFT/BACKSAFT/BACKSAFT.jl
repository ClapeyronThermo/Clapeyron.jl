struct BACKSAFTParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    c::SingleParam{Float64}
    alpha::SingleParam{Float64}
end


abstract type BACKSAFTModel <: SAFTModel end
@newmodel BACKSAFT BACKSAFTModel BACKSAFTParam

export BACKSAFT
function BACKSAFT(components::Array{String,1}; idealmodel=BasicIdeal, userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["SAFT/BACKSAFT","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    segment = params["m"]
    c = params["c"]
    k = params["k"]
    alpha = params["alpha"]
    sigma = params["vol"]
    sigma.values .*= 6/N_A/1e6/π
    sigma.values .^= 1/3
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    packagedparams = BACKSAFTParam(segment, sigma, epsilon, c, alpha)
    references = ["TODO BACKSAFT", "TODO BACKSAFT"]

    return BACKSAFT(packagedparams, idealmodel; references=references, verbose=verbose)
end

include("equations.jl")

#=
function create_BACKSAFTParams(raw_params; combiningrule_ϵ = "Berth")
    like_params_dict, unlike_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "c", "vol", "epsilon", "alpha"];
                     unlike_params = ["k"])

    segment = like_params_dict["m"]

    sigma = like_params_dict["vol"]
    map!(x->(6*x/N_A/1e6/π)^(1/3), values(sigma))
    merge!(sigma, combining_sigma(sigma))

    epsilon = like_params_dict["epsilon"]
    k = unlike_params_dict["k"]
    merge!(epsilon, combining_epsilon(epsilon, sigma, k))

    c = like_params_dict["c"]

    alpha = like_params_dict["alpha"]
    return BACKSAFTParams(segment, sigma, epsilon, c, alpha)
end

include("equations.jl")
=#