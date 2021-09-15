struct sCKSAFTParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type sCKSAFTModel <: CKSAFTModel end
@newmodel sCKSAFT sCKSAFTModel sCKSAFTParam
export sCKSAFT
function sCKSAFT(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, ["SAFT/sCKSAFT","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    segment = params["m"]
    k = params["k"]
    sigma = params["vol"]
    sigma.values .*= 6*0.74048/N_A/1e6/π
    sigma.values .^= 1/3
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]

    packagedparams = sCKSAFTParam(segment, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["TODO sCKSAFT", "TODO sCKSAFT"]

    model = sCKSAFT(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function a_disp(model::sCKSAFTModel, V, T, z)
    ∑z = ∑(z)
    m = model.params.segment.values
    m̄ = dot(z, m)/∑z
    vs = V/(N_A*∑z*m̄)
    it = @comps
    v̄Ȳ = ∑(z[i]*z[j]*m[i]*m[j]*(@f(d,i,j)^3/√(2))*(exp(@f(u,i,j)/T/2)-1) for i ∈ it for j ∈ it)/∑(z[i]*z[j]*m[i]*m[j] for i ∈ it for j ∈ it)
    #res = 36*log(vs/(vs+v̄Ȳ)) 
    #this formulation is prone to horrible floating point issues.
    #writing this term in log1p form solves the problem.
    res = -36*log1p(v̄Ȳ/vs)
    return res
end

function d(model::sCKSAFTModel, V, T, z, i)
    ϵ = model.params.epsilon.diagvalues
    σ = model.params.sigma.diagvalues
    res = σ[i] * (1 - 0.333exp(-3ϵ[i]/T))
    return res
end

function u(model::sCKSAFTModel, V, T, z, i, j)
    ϵ0 = model.params.epsilon.values[i,j]
    res = ϵ0*(1-10/T)
    return res
end

function Δ(model::sCKSAFTModel, V, T, z, i, j, a, b)
    ϵ_assoc = model.params.epsilon_assoc.values[i,j][a,b]
    κ = model.params.bondvol.values[i,j][a,b]
    g = @f(g_hsij,i,j)
    res = g*@f(d,i,j)^3*(exp(ϵ_assoc/T)-1)*κ
    return res
end
