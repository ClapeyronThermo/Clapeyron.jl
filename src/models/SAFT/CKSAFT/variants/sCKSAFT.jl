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
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `vol`: Single Parameter (`Float64`) - Segment Volume `[dm³]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

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

function x0_crit_pure(model::sCKSAFTModel,z)
    T = T_scale(model,z)
    lb_v = lb_volume(model,T,z)/sum(z)
    res = (5.0, log10(lb_v/0.3))
    return res
end

function a_disp(model::sCKSAFTModel, V, T, z, _data = @f(data))
    _d, m̄, ζi, Σz = _data
    m = model.params.segment.values
    m̄1 = dot(z, m)/Σz
    vs = V/(N_A*Σz*m̄)
    
    m̃ = zero(Σz+one(eltype(model)))
    v̄Ȳm̃  = zero(T + Σz+one(eltype(model)))
    sqrt2⁻¹ = 0.7071067811865476 #1/sqrt(2)
    Tinv2 = 1/(2*T)
    for i in @comps
        mi,zi,di = m[i],z[i],_d[i]
        mizi2 = zi*zi*mi*mi
        m̃ += mizi2
        di3 = di*di*di*sqrt2⁻¹
        v̄Ȳm̃ += mizi2*di3*expm1(@f(u,i,i)*Tinv2)
        for j in 1:(i-1)
            mj,zj,dj = m[j],z[j],_d[j]
            mijzij2 = 2*mi*mj*zi*zj
            m̃ += mijzij2
            dij = 0.5*(di+dj)
            dij3 = dij*dij*dij*sqrt2⁻¹
            v̄Ȳm̃ += mijzij2*dij3*expm1(@f(u,i,j)*Tinv2)
        end
    end
    v̄Ȳ = v̄Ȳm̃/m̃
    #it = @comps
    #vy1 = ∑(z[i]*z[j]*m[i]*m[j]*((0.5*_d[i] + 0.5*_d[j])^3/√(2))*(exp(@f(u,i,j)/T/2)-1) for i ∈ it for j ∈ it)
    #m1 = ∑(z[i]*z[j]*m[i]*m[j] for i ∈ it for j ∈ it)
    #@show m1,m̃
    #@show vy1,v̄Ȳm̃
    #v̄Ȳ = vy1/m1
    #res = 36*log(vs/(vs+v̄Ȳ)) 
    #this formulation is prone to horrible floating point issues.
    #writing this term in log1p form solves the problem.
    res = -36*log1p(v̄Ȳ/vs)
    return res
end

d(model::sCKSAFTModel, V, T, z) = ck_diameter(model, T, z, 0.333)

ck_c(model::sCKSAFTModel) = FillArrays.Fill(-10.0,length(model))

function Δ(model::sCKSAFTModel, V, T, z, i, j, a, b, _data = @f(data))
    _d, m̄, ζi, Σz = _data
    ϵ_assoc = model.params.epsilon_assoc.values[i,j][a,b]
    κ = model.params.bondvol.values[i,j][a,b]
    g = @f(g_hs,i,j,_data)
    di,dj = _d[i],_d[j]
    dij = 0.5*(di+dj)
    res = g*dij^3*(exp(ϵ_assoc/T)-1)*κ
    return res
end
