abstract type sPCSAFTModel <: PCSAFTModel end

struct sPCSAFT{T <: IdealModel} <: sPCSAFTModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    sites::SiteParam
    params::PCSAFTParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
    water::SpecialComp
end

@registermodel sPCSAFT

export sPCSAFT
function sPCSAFT(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
    
    params,sites = getparams(components, ["SAFT/PCSAFT/sPCSAFT","properties/molarmass.csv"]; 
    userlocations=userlocations, 
    verbose=verbose,
    ignore_missing_singleparams = ["kT"])
    
    water = SpecialComp(components,["water08"])
    icomponents = 1:length(components)
    segment = params["m"]
    k0 = params["k"]
    n = length(components)
    k1 = get(params,"kT",PairParam("kT",components,zeros(n,n)))
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"])
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]

<<<<<<< HEAD
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon,k0, k1, epsilon_assoc, bondvol)
=======
    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
>>>>>>> 4c2f8f9c7800a07f045f1ae4cb45e1d9d808698a
    references = ["10.1021/ie020753p"]

    model = sPCSAFT(components,icomponents,sites,packagedparams,init_idealmodel,assoc_options,references,water)
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
