
abstract type gcPCPSAFTModel <: PCPSAFTModel end

struct HeterogcPCPSAFTParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    comp_segment::SingleParam{T}
    comp_sigma::PairParam{T}
    comp_epsilon::PairParam{T}
    dipole::SingleParam{T}
    dipole2::SingleParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
end

function HeterogcPCPSAFTParam(Mw,m,σ,ϵ,mc,σc,ϵc,μ,μ2,ϵijab,β)
    return build_parametric_param(HeterogcPCPSAFTParam,Mw,m,σ,ϵ,mc,σc,ϵc,μ,μ2,ϵijab,β)
end

@newmodelgc HeterogcPCPSAFT gcPCPSAFTModel HeterogcPCPSAFTParam{T} true
default_references(::Type{HeterogcPCPSAFT}) = ["10.1021/ie0003887", "10.1021/ie010954d"]
default_locations(::Type{HeterogcPCPSAFT}) = ["SAFT/PCSAFT/gcPCPSAFT/hetero/","properties/molarmass_groups.csv"]
default_gclocations(::Type{HeterogcPCPSAFT}) = ["SAFT/PCSAFT/gcPCPSAFT/hetero/HeterogcPCPSAFT_groups.csv","SAFT/PCSAFT/gcPCPSAFT/hetero/HeterogcPCPSAFT_intragroups.csv"]
default_ignore_missing_singleparams(::Type{HeterogcPCPSAFT}) = String["k"]

function transform_params(::Type{HeterogcPCPSAFT},params,groups)
    components = groups.components
    gc_sigma = params["sigma"]
    gc_sigma.values .*= 1E-10
    gc_epsilon = params["epsilon"]

    #mixing for segment
    gc_segment = params["segment"]
    segment = group_sum(groups,gc_segment)
    params["comp_segment"] = segment

    #mixing for comp_epsilon
    epsilon = group_sum(groups,gc_epsilon .* gc_segment)
    epsilon ./= segment
    epsilon = SingleParam("epsilon",groups.components,epsilon)
    params["comp_epsilon"] = epsilon_LorentzBerthelot(epsilon)

    #mixing for comp_sigma
    gc_sigma = deepcopy(params["sigma"])
    gc_sigma.values .^= 3
    gc_sigma.values .*= gc_segment.values
    sigma = group_sum(groups,gc_sigma)
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    params["comp_sigma"] = sigma_LorentzBerthelot(sigma)

    #mix gc_sigma, gc_epsilon
    params = saft_lorentz_berthelot(params)

    #mix sites
    sites = params["sites"]
    comp_sites = gc_to_comp_sites(sites,groups)
    params["sites"] = comp_sites
    gc_epsilon_assoc = params["epsilon_assoc"]
    gc_bondvol = params["bondvol"]
    assoc_options = params["assoc_options"]
    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,gc_sigma,assoc_options,sites) #combining rules for association
    params["bondvol"] = gc_to_comp_sites(gc_bondvol,comp_sites)
    params["epsilon_assoc"] = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)
    
    #mixing for dipole
    gc_μ = get!(params,"dipole") do
        SingleParam("dipole",components)
    end
    
    gc_μ2 = SingleParam("Dipole squared",groups.flattenedgroups, gc_μ.^2 ./ gc_segment ./ k_B*1e-36*(1e-10*1e-3))
    dipole2 = group_sum(groups,gc_μ2)
    dipole2 = SingleParam("Dipole squared",components, dipole2 ./ segment)
    dipole = SingleParam("Dipole",components, sqrt.(dipole2 .* k_B ./ 1e-36 ./ (1e-10*1e-3)))
    params["dipole"] = dipole
    params["dipole2"] = dipole2
    return params
end

""" 
    gcPCPSAFTModel <: PCSAFTModel

    HeterogcPCPSAFT(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `m`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `k`: Pair Parameter (`Float64`) - Binary Interaction Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models

- `idealmodel`: Ideal Model
## Description

Heterosegmented Group-contribution Polar Perturbed-Chain SAFT (Hetero-gc-PCP-SAFT)

## References
1. Gross, J., Spuhl, O., Tumakaka, F. & Sadowski, G. (2003). Modeling Copolymer Systems Using the Perturbed-Chain SAFT Equation of State. Industrial & Engineering Chemistry Research, 42, 1266-1274. [doi:10.1021/ie020509y](https://doi.org/10.1021/ie020509y)
2. Sauer, E., Stavrou, M. & Gross, J. (2014). Comparison between a Homo- and a Heterosegmented Group Contribution Approach Based on the Perturbed-Chain Polar Statistical Associating Fluid Theory Equation of State. Industrial & Engineering Chemistry Research, 53(38), 14854–14864. [doi:10.1021/ie502203w](https://doi.org/10.1021/ie502203w)

## List of available groups
|Name    |Description         |
|--------|--------------------|
|CH3     |Methyl              |
|CH2     |Methylene           |
|CH      |                    |
|C       |                    |
|CH2=    |Terminal alkene     |
|CH=     |                    |
|=C<     |                    |
|C#CH    |Terminal alkyne     |
|cCH2_pen|Cyclic pentane group|
|cCH_pen |                    |
|cCH2_hex|Cyclic hexane group |
|cCH_hex |                    |
|aCH     |Aromatic group      |
|aCH     |                    |
|OH      |Hydroxyl group      |
|NH2     |Amine group         |
"""
HeterogcPCPSAFT

export HeterogcPCPSAFT

function lb_volume(model::gcPCPSAFTModel, z)
    vk  = model.groups.n_flattenedgroups
    m = model.params.segment.values
    σ = model.params.sigma.values
    σ_idx = linearidx(σ)
    m_idx = linearidx(m)
    val = zero(Base.promote_eltype(m,σ,z))
    for i in 1:length(model)
        val_i = zero(eltype(model))
        vi = vk[i]
        for k in 1:length(model.groups.flattenedgroups)
            mk,σk = m[m_idx[k]],σ[σ_idx[k]]
            val_i += vi[k]*mk*σk*σk*σk
        end
        val += z[i]*val_i
    end
    return π/6*N_A*val
end

function data(model::gcPCPSAFTModel,V,T,z)
    ncomponents = length(model)

    _d = @f(d)
    ζ0,ζ1,ζ2,ζ3 = @f(ζ0123,_d)

    mk = model.params.segment.values
    nk = model.groups.n_flattenedgroups
    m = zero(first(z))*zeros(ncomponents)

    for i in @comps
        m[i] = sum(nk[i].*mk)
    end

    m̄ = dot(z, m)/sum(z)
    return (_d,ζ0,ζ1,ζ2,ζ3,m̄)
end

function a_hc(model::gcPCPSAFTModel, V, T, z,_data=@f(data))
    ngroups = length(model.groups.flattenedgroups)
    _d,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    m = model.params.segment.values
    Σz = sum(z)
    c1 = 1/(1-ζ3)
    c2 = 3ζ2/(1-ζ3)^2
    c3 = 2ζ2^2/(1-ζ3)^3
    if !iszero(ζ3)
        a_hs = bmcs_hs(ζ0,ζ1,ζ2,ζ3)
    else
        a_hs = @f(bmcs_hs_zero_v,_d)
    end
    g_hs = zero(V+T+first(z))*zeros(ngroups,ngroups)
    for k ∈ @groups
        dₖ = _d[k]
        g_hs[k,k] = c1 + dₖ/2*c2 + dₖ^2/4*c3
        for l ∈ @groups
            dₗ = _d[l]
            g_hs[k,l] = c1 + dₖ*dₗ/(dₖ+dₗ)*c2 + (dₖ*dₗ/(dₖ+dₗ))^2*c3
        end
    end
    
    res = zero(V+T+first(z))
    n = model.groups.n_intergroups
    for i ∈ @comps
        for k ∈ @groups
            for l ∈ 1:k
                res+=z[i]*n[i][k,l]*log(g_hs[k,l])
            end
        end
    end
    #return  m̄*@f(a_hs) - ∑(z[i]*(m[i]-1)*log(@f(g_hs,i,i)) for i ∈ @comps)/Σz
    return m̄*a_hs - res/Σz
end

function a_disp(model::gcPCPSAFTModel, V, T, z,_data=@f(data))
    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Σz = sum(z)
    I₁ = @f(I,1,_data)
    I₂ = @f(I,2,_data)
    C₁ = @f(C1,_data)
    m2ϵσ3₁,m2ϵσ3₂ = @f(m2ϵσ3)
    return -2*π*N_A*Σz/V*I₁*m2ϵσ3₁ - π*m̄*N_A*Σz/V*C₁*I₂*m2ϵσ3₂
end

function m2ϵσ3(model::gcPCPSAFTModel, V, T, z)
    n = model.groups.n_flattenedgroups
    m = model.params.segment.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    m2ϵσ3₂ = zero(T+first(z))
    m2ϵσ3₁ = m2ϵσ3₂
    @inbounds for i ∈ @comps
        for j ∈ @comps
            for k ∈ @groups
                for l ∈ @groups
                    constant = z[i]*z[j]*n[i][k]*n[j][l]*m[k]*m[l] * σ[k,l]^3
                    exp1 = (ϵ[k,l]/T)
                    exp2 = exp1*exp1
                    m2ϵσ3₁ += constant*exp1
                    m2ϵσ3₂ += constant*exp2
                end
            end
        end
    end
    Σz = sum(z)
    k = (1/Σz)^2
    return k*m2ϵσ3₁,k*m2ϵσ3₂
    #return ∑(z[i]*z[j]*m[i]*m[j] * (ϵ[i,j]*(1)/T)^n * σ[i,j]^3 for i ∈ @comps, j ∈ @comps)/(sum(z)^2)
end

function d(model::gcPCPSAFTModel, V, T, z)
    ϵᵢᵢ = diagvalues(model.params.epsilon)
    σᵢᵢ = diagvalues(model.params.sigma) 
    return σᵢᵢ .* (1 .- 0.12 .* exp.(-3ϵᵢᵢ ./ T))
end

function ζ0123(model::gcPCPSAFTModel, V, T, z,_d = @f(d))
    m = model.params.segment.values
    n = model.groups.n_flattenedgroups
    ζ0 = zero(V+T+first(z))
    ζ1 = ζ0
    ζ2 = ζ0
    ζ3 = ζ0
    for i ∈ @comps
        for k ∈ @groups
            dₖ = _d[k]
            zᵢmₖ = z[i]*n[i][k]*m[k]
            d1 = dₖ
            d2 = d1*d1
            d3 = d2*d1
            ζ0 += zᵢmₖ
            ζ1 += zᵢmₖ*d1
            ζ2 += zᵢmₖ*d2
            ζ3 += zᵢmₖ*d3
        end
    end
    NV = N_A*π/6/V
    ζ0 *= NV
    ζ1 *= NV
    ζ2 *= NV
    ζ3 *= NV
    return ζ0,ζ1,ζ2,ζ3
end

function Δ(model::gcPCPSAFTModel, V, T, z, i, j, a, b,_data=@f(data))
    _0 = zero(V+T+first(z))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    κijab = κ[i,j][a,b] 
    iszero(κijab) && return _0
    σ = model.params.sigma.values
    k,l = get_group_idx(model,i,j,a,b)
    gkl = @f(g_hs,k,l,_data)
    res = gkl*σ[k,l]^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κijab
    return res
end

function  Δ(model::HeterogcPCPSAFT, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    Δout = assoc_similar(κ,typeof(V+T+first(z)))
    Δout.values .= false  #fill with zeros, maybe it is not necessary?
    for (idx,(i,j),(a,b)) in indices(Δout)
        k,l = get_group_idx(model,i,j,a,b)
        gkl = @f(g_hs,k,l,_data)
        Δout[idx] = gkl*σ[k,l]^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]
    end
    return Δout
end

#polar overloads
@inline pcp_sigma(model::HeterogcPCPSAFT) = model.params.comp_sigma.values
@inline pcp_epsilon(model::HeterogcPCPSAFT) = model.params.comp_epsilon.values
@inline pcp_segment(model::HeterogcPCPSAFT) = model.params.comp_segment.values
