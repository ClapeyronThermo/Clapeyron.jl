struct sCKSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type sCKSAFTModel <: CKSAFTModel end
@newmodel sCKSAFT sCKSAFTModel sCKSAFTParam
default_references(::Type{sCKSAFT}) = ["10.1021/IE00107A014", "10.1021/ie00056a050","10.1021/ie00044a042"]
default_locations(::Type{sCKSAFT}) = ["SAFT/sCKSAFT","properties/molarmass.csv"]
function transform_params(::Type{sCKSAFT},params)
    k = get(params,"k",nothing)
    sigma = params["vol"]
    sigma.values .*= 6*0.74048/N_A/1e6/π
    sigma.values .^= 1/3
    epsilon = params["epsilon"]
    params["sigma"] = sigma_LorentzBerthelot(sigma)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon, k)
    return params
end

export sCKSAFT

"""
    sCKSAFTModel <: CKSAFTModel

    sCKSAFT(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `vol`: Single Parameter (`Float64`) - Segment Volume [`dm^3`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model

## Description

Simplified Chen and Kreglewski SAFT (sCK-SAFT)

## References
1. Huang, S. H., & Radosz, M. (1990). Equation of state for small, large, polydisperse, and associating molecules. Industrial & Engineering Chemistry Research, 29(11), 2284–2294. [doi:10.1021/ie00107a014](https://doi.org/10.1021/ie00107a014)
2. Huang, S. H., & Radosz, M. (1991). Equation of state for small, large, polydisperse, and associating molecules: extension to fluid mixtures. Industrial & Engineering Chemistry Research, 30(8), 1994–2005. [doi:10.1021/ie00056a050](https://doi.org/10.1021/ie00056a050)
3. Fu, Y.-H., & Sandler, S. I. (1995). A simplified SAFT equation of state for associating compounds and mixtures. Industrial & Engineering Chemistry Research, 34(5), 1897–1909. [doi:10.1021/ie00044a042](https://doi.org/10.1021/ie00044a042)
"""
sCKSAFT

function x0_crit_pure(model::sCKSAFTModel)
    lb_v = lb_volume(model)
    res = (5.0, log10(lb_v/0.3))
    return res
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
    ϵ = model.params.epsilon.values[i,i]
    σ = model.params.sigma.values[i,i]
    res = σ * (1 - 0.333exp(-3ϵ/T))
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
