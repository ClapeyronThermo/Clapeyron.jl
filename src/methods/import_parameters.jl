function create_PCSAFTParams(raw_params)
    pure_params_dict, pair_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon", "n_H", "n_e"];
                     pair_params = ["k"], assoc_params = ["epsilon_assoc", "bond_vol"])
    segment = pure_params_dict["m"]
    sigma = pure_params_dict["sigma"] # If multiplication by a factor is needed, do it here
    # We can also add binary parameters using push!(sigma, Set([x1, x2]) => value) later
    # using a combining rule.
    epsilon = pure_params_dict["epsilon"]
    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    k = pair_params_dict["k"]
    return PCSAFTParams(segment, sigma, epsilon, epsilon_assoc, bond_vol, k)
end
    
function create_SAFTVRMieParams(raw_params)
    pure_params_dict, pair_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon", "lambdaA", "lambdaR"])
    segment = pure_params_dict["m"]
    sigma = pure_params_dict["sigma"]
    epsilon = pure_params_dict["epsilon"]
    lambdaA = pure_params_dict["lambdaA"]
    lambdaR = pure_params_dict["lambdaR"]
    return SAFTVRMieParams(segment, sigma, epsilon, lambdaA, lambdaR)
end
    
function create_ogSAFTParams(raw_params)
    pure_params_dict, pair_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon"])
    segment = pure_params_dict["m"]
    sigma = pure_params_dict["sigma"]
    epsilon = pure_params_dict["epsilon"]
    return ogSAFTParams(segment, sigma, epsilon)
end
