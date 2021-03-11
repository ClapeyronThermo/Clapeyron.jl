struct CKSAFTParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    c::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end


abstract type CKSAFTModel <: SAFTModel end
@newmodel CKSAFT CKSAFTModel CKSAFTParam

export CKSAFT
function CKSAFT(components::Array{String,1}; idealmodel=BasicIdeal, userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["SAFT/CKSAFT","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    segment = params["m"]
    c = params["c"]
    k = params["k"]
    sigma = params["vol"]
    sigma.values .*= 6*0.74048/N_A/1e6/π
    sigma.values .^= 1/3
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bond_vol"]
    sites = SiteParam(Dict("e" => params["n_e"], "H" => params["n_H"]))

    packagedparams = CKSAFTParam(segment, sigma, epsilon,c, epsilon_assoc, bondvol)
    references = ["TODO CKSAFT", "TODO CKSAFT"]

    return CKSAFT(packagedparams, sites, idealmodel; references=references, verbose=verbose)
end

include("equations.jl")

