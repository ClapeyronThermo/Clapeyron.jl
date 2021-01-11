struct SAFTgammaMieParam <: EoSParam
    segment::SingleParam{Int}
    shapefactor::SingleParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type SAFTgammaMieModel <: SAFTModel end
@newmodelgc SAFTgammaMie SAFTgammaMieModel SAFTgammaMieParam

export SAFTgammaMie
function SAFTgammaMie(components::Array{<:Any,1}; idealmodels::Array{String,1}=String[], userlocations::Array{String,1}=String[], verbose=false)
    groups = buildspecies(components, ["SAFT/SAFTgammaMie/groups.csv"]; verbose=verbose)
    params = getparams(groups, ["SAFT/SAFTgammaMie"]; userlocations=userlocations, modelname="SAFTgammaMie", verbose=verbose)

    segment = params["vst"]
    shapefactor = params["S"]

    params["sigma"].values .*= 1E-10
    sigma = combining_sigma(params["sigma"])
    epsilon = combining_epsilon(params["epsilon"], sigma; rules="Hudson-McCoubrey")
    lambda_a = combining_lambda(params["lambda_a"])
    lambda_r = combining_lambda(params["lambda_r"])

    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bond_vol"]

    sites = getsites(Dict("e1" => params["n_e1"], "e2" => params["n_e2"], "H" => params["n_H"]))

    packagedparams = SAFTgammaMieParam(segment, shapefactor, lambda_a, lambda_r, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1063/1.4851455", "10.1021/je500248h"]

    return SAFTgammaMie(packagedparams, groups, sites; references=references)
end

include("equations.jl")
