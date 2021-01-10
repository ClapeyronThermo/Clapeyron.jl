struct PCSAFTParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type PCSAFTModel <: SAFTModel end
@newmodel PCSAFT PCSAFTModel PCSAFTParam

export PCSAFT
function PCSAFT(components::Array{String,1}; idealmodels::Array{String,1}=String[], userlocations::Array{String,1}=String[], combiningrule_epsilon="Berth", verbose=false)
    params = getparams(components, ["SAFT/PCSAFT"]; userlocations=userlocations, modelname="PCSAFT", verbose=verbose)
    # ideal = getideal(components, idealmodels)
    segment = params["m"]
    k = params["k"]
    params["sigma"].values .*= 1E-10
    sigma = combining_sigma(params["sigma"])
    epsilon = combining_epsilon(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bond_vol"]
    sites = getsites(Dict("e" => params["n_e"], "H" => params["n_H"]))

    packagedparams = PCSAFTParam(segment, sigma, epsilon, epsilon_assoc, bondvol)

    return PCSAFT(packagedparams, sites)
end

include("equations.jl")
