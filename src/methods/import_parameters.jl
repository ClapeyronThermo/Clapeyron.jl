function create_PCSAFTParams(raw_params; combiningrule_ϵ = "Berth")
    pure_params_dict, pair_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon", "n_H", "n_e"];
                     pair_params = ["k"], assoc_params = ["epsilon_assoc", "bond_vol"])
    segment = pure_params_dict["m"]
    sigma = pure_params_dict["sigma"]
    map!(x->x*1E-10, values(sigma))
    # We can also add binary parameters using push!(sigma, Set([x1, x2]) => value) later
    # using a combining rule.
    epsilon = pure_params_dict["epsilon"]
    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(pure_params_dict["n_e"])
        n_sites[i] = Dict()
        n_sites[i][Set(["e"])] = pure_params_dict["n_e"][i]
        n_sites[i][Set(["H"])] = pure_params_dict["n_H"][i]
    end
    k = pair_params_dict["k"]
    return PCSAFTParams(segment, sigma, epsilon, epsilon_assoc, bond_vol, n_sites, k)
end

function create_SAFTVRMieParams(raw_params)
    pure_params_dict, pair_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon", "lambdaA", "lambdaR","n_H","n_e"];
                     assoc_params = ["epsilon_assoc", "bond_vol"])
    segment = pure_params_dict["m"]
    sigma = pure_params_dict["sigma"]
    map!(x->x*1E-10, values(sigma))
    epsilon = pure_params_dict["epsilon"]
    lambdaA = pure_params_dict["lambdaA"]
    lambdaR = pure_params_dict["lambdaR"]
    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(pure_params_dict["n_H"])
        n_sites[i] = Dict()
        n_sites[i][Set(["e"])] = pure_params_dict["n_e"][i]
        n_sites[i][Set(["H"])] = pure_params_dict["n_H"][i]
    end
    return SAFTVRMieParams(segment, sigma, epsilon, lambdaA, lambdaR, epsilon_assoc, bond_vol, n_sites)
end

function create_ogSAFTParams(raw_params)
    pure_params_dict, pair_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon","n_H","n_e"];
                     assoc_params = ["epsilon_assoc", "bond_vol"])
    segment = pure_params_dict["m"]
    sigma = pure_params_dict["sigma"]
    map!(x->x*1E-10, values(sigma))
    epsilon = pure_params_dict["epsilon"]
    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(pure_params_dict["n_H"])
        n_sites[i] = Dict()
        n_sites[i][Set(["e"])] = pure_params_dict["n_e"][i]
        n_sites[i][Set(["H"])] = pure_params_dict["n_H"][i]
    end
    return ogSAFTParams(segment, sigma, epsilon, epsilon_assoc, bond_vol, n_sites)
end
