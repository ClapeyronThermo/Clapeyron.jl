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
        n_sites[i]["e"] = like_params_dict["n_e"][i]
        n_sites[i]["H"] = like_params_dict["n_H"][i]
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
        n_sites[i]["e"] = like_params_dict["n_e"][i]
        n_sites[i]["H"] = like_params_dict["n_H"][i]
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
        n_sites[i]["e"] = like_params_dict["n_e"][i]
        n_sites[i]["H"] = like_params_dict["n_H"][i]
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
        n_sites[i]["e"] = like_params_dict["n_e"][i]
        n_sites[i]["H"] = like_params_dict["n_H"][i]
    end
    return ogSAFTParams(segment, sigma, epsilon, epsilon_assoc, bond_vol, n_sites)
end

function create_SAFTgammaMie(raw_params)
    like_params_dict, unlike_params_dict, assoc_params_dict =
        filterparams(raw_params, ["vst", "S", "sigma", "epsilon", "lambda_a", "lambda_r","n_H","n_e1", "n_e2"];
                     unlike_params=["epsilon","lambda_r"], assoc_params = ["epsilon_assoc", "bond_vol"])
    segment = like_params_dict["vst"]
    shapefactor = like_params_dict["S"]

    sigma = like_params_dict["sigma"]
    merge!(sigma, combining_sigma(sigma))
    map!(x->x*1E-10, values(sigma))

    epsilon = like_params_dict["epsilon"]
    merge!(epsilon, unlike_params_dict["epsilon"])
    # may need modifications
    merge!(epsilon, combining_epsilon(epsilon, sigma, Dict();rules_no_k = "Hudson-McCoubrey"))

    lambda_a = like_params_dict["lambda_a"]
    merge!(lambda_a, combining_lambda(lambda_a))

    lambda_r = like_params_dict["lambda_r"]
    merge!(lambda_r, unlike_params_dict["lambda_r"])
    merge!(lambda_r, combining_lambda(lambda_r))

    epsilon_assoc = assoc_params_dict["epsilon_assoc"]
    bond_vol = assoc_params_dict["bond_vol"]
    n_sites = Dict()
    for i in keys(like_params_dict["n_H"])
        n_sites[i] = Dict()
        n_sites[i]["e1"] = like_params_dict["n_e1"][i]
        n_sites[i]["e2"] = like_params_dict["n_e2"][i]
        n_sites[i]["H"] = like_params_dict["n_H"][i]
    end
    return SAFTgammaMieParams(segment, shapefactor, lambda_a, lambda_r, sigma, epsilon, epsilon_assoc, bond_vol, n_sites)
end

function create_vdWParams(raw_params)
    like_params_dict, unlike_params_dict, assoc_params_dict =
        filterparams(raw_params, ["a_vdW", "b_vdW"];
                     unlike_params = ["k"])

    b_vdW = like_params_dict["b_vdW"]

    merge!(b_vdW, combining_sigma(b_vdW))

    a_vdW = like_params_dict["a_vdW"]
    k = unlike_params_dict["k"]
    merge!(epsilon, combining_epsilon(a_vdW, b_vdW, k))

    return vdWParams(a_vdW,b_vdW)
end
