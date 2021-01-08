# include("equations")

struct PCSAFTParam <: EoSParam
    segment::SingleParam
    sigma::PairParam
    epsilon::PairParam
    epsilon_assoc::AssocParam
    bondvol::AssocParam
end

struct PCSAFTModel <: SAFTModel
    components::Array{String,1}
    sites::Array{Array{String,1},1}
    nsites::Array{Array{Int,1},1}
    params::PCSAFTParam
end

export PCSAFT
function PCSAFT(components::Array{String,1}; idealmodels::Array{String,1}=String[], usermodels::Array{String,1}=String[], combiningrule_epsilon="Berth", verbose=false)
    params, sites = getparams(components, ["SAFT/PCSAFT"]; usermodels=usermodels, modelname="PCSAFT", verbose=verbose)
    # ideal = getideal(components, idealmodels)
    segment = params["m"]
    k = params["k"]
    params["sigma"].values .*= 1E-10
    sigma = combining_sigma(params["sigma"])
    epsilon = combining_epsilon(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bond_vol"]
    nsites = getnsites(sites, Dict("e" => params["n_e"], "H" => params["n_H"]))
    return PCSAFTModel(components, sites, nsites, PCSAFTParam(segment, sigma, epsilon, epsilon_assoc, bondvol))
end

function getnsites(sites::Array{Array{String,1},1}, pairs::Dict{String,SingleParam{Int}})
    numberofcomponents = length(sites)
    nsites = [Int[] for x in 1:numberofcomponents]
    for i ∈ 1:numberofcomponents, j ∈ 1:length(sites[i])
        append!(nsites[i], pairs[sites[i][j]].values[i])
    end
    return nsites
end


