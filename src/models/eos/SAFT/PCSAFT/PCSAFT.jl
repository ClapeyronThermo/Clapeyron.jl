struct PCSAFTParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type PCSAFTModel <: NonGCSAFTModel end
@newmodel PCSAFT PCSAFTModel PCSAFTParam

export PCSAFT
function PCSAFT(components::Array{String,1}; idealmodel::String="", userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["SAFT/PCSAFT"]; userlocations=userlocations, verbose=verbose)
    segment = params["m"]
    k = params["k"]
    params["sigma"].values .*= 1E-10
    sigma = combining_sigma(params["sigma"])
    epsilon = combining_epsilon(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    sites = getsites(Dict("e" => params["n_e"], "H" => params["n_H"]))

    packagedparams = PCSAFTParam(segment, sigma, epsilon, epsilon_assoc, bondvol)
    idealmodel = idealmodelselector(idealmodel, components)
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]

    return PCSAFT(packagedparams, sites, idealmodel; references=references)
end

include("equations.jl")
