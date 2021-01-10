abstract type sPCSAFTModel <: PCSAFTModel end
@newmodel sPCSAFT sPCSAFTModel PCSAFTParam

export sPCSAFT
function sPCSAFT(components::Array{String,1}; idealmodels::Array{String,1}=String[], userlocations::Array{String,1}=String[], combiningrule_epsilon="Berth", verbose=false)
    params = getparams(components, ["SAFT/PCSAFT", "SAFT/sPCSAFT"]; userlocations=userlocations, modelname="sPCSAFT", verbose=verbose)
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

    return sPCSAFT(packagedparams, sites)
end

function a_hc(model::sPCSAFTModel, z, V, T)
    x = z/sum(z)
    m = model.params.segment.values
    m̄ = ∑(x .* m)
    return m̄*@f(a_hs) - (m̄-1)*log(@f(g_hs))
end

function g_hs(model::sPCSAFTModel, z, V, T)
    η = @f(ζ,3)
    return (1-η/2)/(1-η)^3
end

function a_hs(model::sPCSAFTModel, z, V, T)
    η = @f(ζ,3)
    return (4η-3η^2)/(1-η)^2
end

function Δ(model::sPCSAFTModel, z, V, T, i, j, a, b)
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    g_hs_ = @f(g_hs)
    return g_hs_*σij^3*(exp(ϵ_associjab/T)-1)*κijab
end
