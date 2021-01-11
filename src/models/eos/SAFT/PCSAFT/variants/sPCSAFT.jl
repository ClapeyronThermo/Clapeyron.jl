abstract type sPCSAFTModel <: PCSAFTModel end
@newmodel sPCSAFT sPCSAFTModel PCSAFTParam

export sPCSAFT
function sPCSAFT(components::Array{String,1}; idealmodel::String="", userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["SAFT/PCSAFT", "SAFT/PCSAFT/sPCSAFT"]; userlocations=userlocations, modelname="sPCSAFT", verbose=verbose)
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
    references = ["10.1021/ie020753p"]

    return sPCSAFT(packagedparams, sites, idealmodel; references=references)
end

function a_hc(model::sPCSAFTModel, V, T, z)
    x = z/sum(z)
    m = model.params.segment.values
    m̄ = ∑(x .* m)
    return m̄*@f(a_hs) - (m̄-1)*log(@f(g_hs))
end

function g_hs(model::sPCSAFTModel, V, T, z)
    η = @f(ζ,3)
    return (1-η/2)/(1-η)^3
end

function a_hs(model::sPCSAFTModel, V, T, z)
    η = @f(ζ,3)
    return (4η-3η^2)/(1-η)^2
end

function Δ(model::sPCSAFTModel, V, T, z, i, j, a, b)
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    g_hs_ = @f(g_hs)
    return g_hs_*σij^3*(exp(ϵ_associjab/T)-1)*κijab
end
