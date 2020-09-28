function create_PCSAFTParams(raw_params; combiningrule_ϵ = "Berth")
    like_params_dict, unlike_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon", "n_H", "n_e"];
                     unlike_params = ["k"], assoc_params = ["epsilon_assoc", "bond_vol"])

    segment = like_params_dict["m"]

    sigma = like_params_dict["sigma"]
    map!(x->x*1E-10, values(sigma))
    merge!(sigma, combining_sigma(sigma))

    epsilon = like_params_dict["epsilon"]
    k = unlike_params_dict["k"]
    merge!(epsilon, combining_epsilon(epsilon, sigma, k))

    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(like_params_dict["n_e"])
        n_sites[i] = Dict()
        n_sites[i][Set(["e"])] = like_params_dict["n_e"][i]
        n_sites[i][Set(["H"])] = like_params_dict["n_H"][i]
    end

    return PCSAFTParams(segment, sigma, epsilon, epsilon_assoc, bond_vol, n_sites)
end

function create_sPCSAFTParams(raw_params; combiningrule_ϵ = "Berth")
    like_params_dict, unlike_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon", "n_H", "n_e"];
                     unlike_params = ["k"], assoc_params = ["epsilon_assoc", "bond_vol"])

    segment = like_params_dict["m"]

    sigma = like_params_dict["sigma"]
    map!(x->x*1E-10, values(sigma))
    merge!(sigma, combining_sigma(sigma))

    epsilon = like_params_dict["epsilon"]
    k = unlike_params_dict["k"]
    merge!(epsilon, combining_epsilon(epsilon, sigma, k))

    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(like_params_dict["n_e"])
        n_sites[i] = Dict()
        n_sites[i][Set(["e"])] = like_params_dict["n_e"][i]
        n_sites[i][Set(["H"])] = like_params_dict["n_H"][i]
    end

    return sPCSAFTParams(segment, sigma, epsilon, epsilon_assoc, bond_vol, n_sites, k)
end

function create_SAFTVRMieParams(raw_params)
    like_params_dict, unlike_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon", "lambdaA", "lambdaR","n_H","n_e"];
                     unlike_params=["epsilon","lambdaR"],assoc_params = ["epsilon_assoc", "bond_vol"])
    segment = like_params_dict["m"]

    sigma = like_params_dict["sigma"]
    merge!(sigma, combining_sigma(sigma))
    map!(x->x*1E-10, values(sigma))

    epsilon = like_params_dict["epsilon"]
    merge!(epsilon, unlike_params_dict["epsilon"])
    merge!(epsilon, combining_epsilon(epsilon, sigma, Dict();rules_no_k = "Hudson-McCoubrey"))

    lambdaA = like_params_dict["lambdaA"]
    merge!(lambdaA, combining_lambda(lambdaA))

    lambdaR = like_params_dict["lambdaR"]
    merge!(lambdaR, unlike_params_dict["lambdaR"])
    merge!(lambdaR, combining_lambda(lambdaR))

    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(like_params_dict["n_H"])
        n_sites[i] = Dict()
        n_sites[i][Set(["e"])] = like_params_dict["n_e"][i]
        n_sites[i][Set(["H"])] = like_params_dict["n_H"][i]
    end
    return SAFTVRMieParams(segment, sigma, epsilon, lambdaA, lambdaR, epsilon_assoc, bond_vol, n_sites)
end

function create_ogSAFTParams(raw_params)
    like_params_dict, unlike_params_dict, assoc_params_dict =
        filterparams(raw_params, ["m", "sigma", "epsilon","n_H","n_e"];
                     assoc_params = ["epsilon_assoc", "bond_vol"])
    segment = like_params_dict["m"]
    sigma = like_params_dict["sigma"]
    map!(x->x*1E-10, values(sigma))
    merge!(sigma, combining_sigma(sigma))
    epsilon = like_params_dict["epsilon"]
    merge!(epsilon, combining_epsilon(epsilon, sigma, Dict()))
    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(like_params_dict["n_H"])
        n_sites[i] = Dict()
        n_sites[i][Set(["e"])] = like_params_dict["n_e"][i]
        n_sites[i][Set(["H"])] = like_params_dict["n_H"][i]
    end
    return ogSAFTParams(segment, sigma, epsilon, epsilon_assoc, bond_vol, n_sites)
end
