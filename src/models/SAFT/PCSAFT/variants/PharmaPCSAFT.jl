struct pharmaPCSAFTParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    k::PairParam{T}
    kT::PairParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
    water::SpecialComp
end

abstract type pharmaPCSAFTModel <: PCSAFTModel end
@newmodel pharmaPCSAFT pharmaPCSAFTModel pharmaPCSAFTParam{T}

function pharmaPCSAFTParam(Mw,m,σ,ϵ,k,kT,ϵijab,β,water)
    return build_parametric_param(pharmaPCSAFTParam,Mw,m,σ,ϵ,k,kT,ϵijab,β,water)
end

default_references(::Type{pharmaPCSAFT}) = ["10.1021/ie0003887", "10.1021/ie010954d","10.1016/j.cep.2007.02.034"]
default_locations(::Type{pharmaPCSAFT}) = ["SAFT/PCSAFT","SAFT/PCSAFT/pharmaPCSAFT","properties/molarmass.csv"]
default_assoc_options(::Type{pharmaPCSAFT}) = AssocOptions(combining = :elliott_runtime)
function transform_params(::Type{pharmaPCSAFT},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    params["kT"] = get(params,"kT",PairParam("kT",components,zeros(length(components))))
    water = SpecialComp(components,["water08","water"])
    params["water"] = water

    if water[] != 0
        assoc = params["assoc_options"]
        if assoc.combining ∈ (:esd,:elliott)
            throw(ArgumentError("pharmaPCSAFT does not support $(assoc.combining) rule with a temperature-dependent sigma. try using :esd_runtime instead"))
        end
    end
    #k needs to be fully instantiated here, because it is stored as a parameter
    params["k"] = get(params,"k",PairParam("k",components,zeros(length(components))))
    sigma,epsilon = params["sigma"],params["epsilon"]
    params["sigma"] = sigma_LorentzBerthelot(sigma)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon) #the mixing is done at runtime
    return params
end

export pharmaPCSAFT

"""
    pharmaPCSAFTModel <: PCSAFTModel

    pharmaPCSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Constant binary Interaction Parameter (no units)
- `kT`: Pair Parameter (`Float64`) - T-dependent Binary Interaction Parameter `[K⁻¹]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`
## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Constant binary Interaction Parameter (no units)
- `kT`: Pair Parameter (`Float64`) - T-dependent Binary Interaction Parameter `[K⁻¹]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`
## Input models
- `idealmodel`: Ideal Model
## Description
Perturbed-Chain SAFT (PC-SAFT), with T dependent kij and water correlation [2] for segment diameter.
For using the water's sigma correlation, `water08` should be selected instead of `water`.
## References
1. Paus, R., Ji, Y., Vahle, L., & Sadowski, G. (2015). Predicting the solubility advantage of amorphous pharmaceuticals: A novel thermodynamic approach. Molecular Pharmaceutics, 12(8), 2823–2833. [doi:10.1021/mp500824d](https://doi.org/10.1021/mp500824d)
2. Cameretti, L. F., & Sadowski, G. (2008). Modeling of aqueous amino acid and polypeptide solutions with PC-SAFT. Genie Des Procedes [Chemical Engineering and Processing], 47(6), 1018–1025. [doi:10.1016/j.cep.2007.02.034](https://doi.org/10.1016/j.cep.2007.02.034)
"""
pharmaPCSAFT

#Δσh20(T) = σ[T] - σconstant
# σ = σconstant + Δσh20(T)
#https://doi.org/10.1016/j.cep.2007.02.034
Δσh20(T) = (10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10
@inline water08_k(model::PCSAFTModel) = 0
@inline water08_k(model::pharmaPCSAFTModel) = model.params.water[]

x0_volume_liquid(model::PCSAFTModel,p,T,z) = lb_volume(model,T,z)*1.7

lb_volume(model::pharmaPCSAFTModel, z) = lb_volume(model,298.15,z)

function lb_volume(model::pharmaPCSAFTModel,T, z)
    seg = model.params.segment.values
    σ = model.params.sigma.values
    k = water08_k(model)
    if k > 0
        Δσ = Δσh20(T)
    else
        Δσ = zero(T)
    end
    val = zero(Base.promote_eltype(model,T,z))
    for i in @comps
        σi = σ[i,i] + (k==i)*Δσ
        val += z[i]*seg[i]*σi*σi*σi
    end
    lb_v = π/6*N_A*val
    return lb_v
end

function d(model::pharmaPCSAFTModel, V, T, z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    _d = zeros(typeof(T*one(eltype(model))),length(z))
    Δσ = Δσh20(T)
    k = water08_k(model)
    for i ∈ @comps
        σᵢ = σ[i,i] + (k==i)*Δσ
        _d[i] = σᵢ*(1 - 0.12*exp(-3ϵ[i,i]/T))
    end

    return _d
end

function d(model::pharmaPCSAFT, V, T, z::SingleComp)
    ϵ = only(model.params.epsilon.values)
    σ = only(model.params.sigma.values)
    k = water08_k(model)
    if k == 1
        σ += Δσh20(T)
    end
    return SA[σ*(1 - 0.12*exp(-3ϵ/T))]
end

function m2ϵσ3(model::pharmaPCSAFTModel, V, T, z)
    m = model.params.segment.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    k0 = model.params.k.values
    k1 = model.params.kT.values
    m2ϵσ3₂ = zero(V+T+first(z))
    m2ϵσ3₁ = m2ϵσ3₂

    k = water08_k(model)
    Δσ = Δσh20(T)
    @inbounds for i ∈ @comps
        zi = z[i]
        mi = m[i]
        ki = (k != i)
        σii = σ[i,i] + (k==i)*Δσ
        constant = zi*zi*mi*mi*(σii)^3
        exp1 = (ϵ[i,i]/T)
        exp2 = exp1*exp1
        m2ϵσ3₁ += constant*exp1
        m2ϵσ3₂ += constant*exp2
        for j ∈ 1:(i-1)
            σij = σ[i,j] + (0.5*(k==i)  +  0.5*(k==j))*Δσ
            constant = zi*z[j]*mi*m[j] * σij^3
            exp1 = ϵ[i,j]*(1 - k0[i,j] - k1[i,j]*T)/T
            exp2 = exp1*exp1
            m2ϵσ3₁ += 2*constant*exp1
            m2ϵσ3₂ += 2*constant*exp2
        end
    end
    Σz = sum(z)
    invn2 = (1/Σz)^2
    return invn2*m2ϵσ3₁,invn2*m2ϵσ3₂
    #return ∑(z[i]*z[j]*m[i]*m[j] * (ϵ[i,j]*(1)/T)^n * σ[i,j]^3 for i ∈ @comps, j ∈ @comps)/(sum(z)^2)
end

function Δ(model::pharmaPCSAFTModel, V, T, z, i, j, a, b,_data=@f(data))
    _0 = zero(V+T+first(z))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    κijab = κ[i,j][a,b]
    iszero(κijab) && return _0
    σ = model.params.sigma.values
    gij = @f(g_hs,i,j,_data)
    k = water08_k(model)
    Δσ = Δσh20(T)
    σij = σ[i,j] + (0.5*(k==i) + 0.5*(k==j))*Δσ
    res = gij*σij^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κijab
    return res
end

function  Δ(model::pharmaPCSAFTModel, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    k = water08_k(model)
    k = model.params.water[]
    Δσ = Δσh20(T)
    Δout = assoc_similar(κ,typeof(V+T+first(z)))
    Δout.values .= false #fill with zeros, maybe it is not necessary?
    for (idx,(i,j),(a,b)) in indices(Δout)
        κijab = κ[idx]
        iszero(κijab) && continue
        gij = @f(g_hs,i,j,_data)
        σij = σ[i,j] + (0.5*(k==i)  +  0.5*(k==j))*Δσ
        Δout[idx] = gij*σij^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κijab
    end
    return Δout
end