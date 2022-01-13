abstract type sPCSAFTModel <: PCSAFTModel end
@newmodel sPCSAFT sPCSAFTModel PCSAFTParam

export sPCSAFT
function sPCSAFT(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
    
    params,sites = getparams(components, ["SAFT/PCSAFT", "SAFT/PCSAFT/sPCSAFT"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["m"]
    k = params["k"]
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1021/ie020753p"]

    model = sPCSAFT(packagedparams, sites, idealmodel; ideal_userlocations, references, verbose, assoc_options)
    return model
end

function a_hc(model::sPCSAFTModel, V, T, z , _data = @f(data))
    _,_,_,_,η,m̄ = _data
    g_hs = (1-η/2)/(1-η)^3
    a_hs = (4η-3η^2)/(1-η)^2
    return m̄*a_hs - (m̄-1)*log(g_hs)
end

function g_hs(model::sPCSAFTModel, V, T, z,_data = @f(data))
    _,_,_,_,η,_ = _data
    return (1-η/2)/(1-η)^3
end

function a_hs(model::sPCSAFTModel, V, T, z)
    _,_,_,_,η,_ = _data
    return (4η-3η^2)/(1-η)^2
end

function Δ(model::sPCSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    g_hs_ = @f(g_hs,_data)
    return g_hs_*σij^3*(exp(ϵ_associjab/T)-1)*κijab
end
