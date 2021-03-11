
abstract type sCKSAFTModel <: CKSAFTModel end


@newmodel sCKSAFT CKSAFTModel CKSAFTParam

export sCKSAFT
function sCKSAFT(components::Array{String,1}; idealmodel=BasicIdeal, userlocations::Array{String,1}=String[], verbose=false)
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
    references = ["TODO sCKSAFT", "sTODO CKSAFT"]

    return sCKSAFT(packagedparams, sites, idealmodel; references=references, verbose=verbose)
end

function a_disp(model::sCKSAFTFamily, V, T, z)
    ∑z = ∑(z)
    x = z.*(1/∑z)
    m = model.params.segment.values
    m̄ = ∑(x[i]*m[i] for i in @comps)
    vs = V/(N_A*∑z*m̄)
    it = @comps
    v̄Ȳ = sum(x[i]*x[j]*m[i]*m[j]*(@f(d,i,j)^3/√(2))*(exp(@f(u,i,j)/T/2)-1) for i in it for j in it)/sum(x[i]*x[j]*m[i]*m[j] for i in it for j in it)
    return 36*log(vs/(vs+v̄Ȳ))
end

function d(model::sCKSAFTFamily, V, T, z, component)
    ϵ = model.params.epsilon.values[component]
    σ = model.params.sigma.values[component]
    return σ * (1 - 0.333exp(-3ϵ/T))
end

function u(model::sCKSAFTFamily, V, T, z, i,j)
    ϵ0 = model.params.epsilon.values[i,j]
    return ϵ0*(1-10/T)
end

function Δ(model::sCKSAFTFamily, V, T, z, i, j, a, b)
    ϵ_assoc = model.params.epsilon_assoc.values[i,j][a,b]
    κ = model.params.bondvol.values[i,j][a,b]
    g = @f(g_hsij,i,j)
    return g*d(model,z,v,T,i,j)^3*(exp(ϵ_assoc/T)-1)*κ
end
