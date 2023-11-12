struct gcPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type gcPCSAFTModel <: PCSAFTModel end
@newmodelgc gcPCSAFT gcPCSAFTModel gcPCSAFTParam true true
default_references(::Type{gcPCSAFT}) = ["10.1021/ie0003887", "10.1021/ie010954d"]
default_locations(::Type{gcPCSAFT}) = ["SAFT/PCSAFT/gcPCSAFT","properties/molarmass_groups.csv"]
default_gclocations(::Type{gcPCSAFT}) = ["SAFT/PCSAFT/gcPCSAFT/gcPCSAFT_groups.csv","SAFT/PCSAFT/gcPCSAFT/gcPCSAFT_intragroups.csv"]
function transform_params(::Type{gcPCSAFT},params,groups)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    params = saft_lorentz_berthelot(params)
    
    sites = params["sites"]
    comp_sites = gc_to_comp_sites(sites,groups)
    params["sites"] = comp_sites
    gc_epsilon_assoc = params["epsilon_assoc"]
    gc_bondvol = params["bondvol"]
    assoc_options = params["assoc_options"]
    
    sigma = params["sigma"]
    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,sigma,assoc_options) #combining rules for association
    params["bondvol"] = gc_to_comp_sites(gc_bondvol,comp_sites)
    params["epsilon_assoc"] = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)
    return params
end

"""
    gcPCSAFTModel <: PCSAFTModel
    gcPCSAFT(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
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
Heterogeneous Group-contribution Perturbed-Chain SAFT (gc-PC-SAFT)
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
gcPCSAFT

export gcPCSAFT

function lb_volume(model::gcPCSAFTModel, z = SA[1.0])
    vk  = model.groups.n_flattenedgroups
    seg = model.params.segment.values
    σ = model.params.sigma.values
    val = π/6*N_A*sum(z[i]*sum(vk[i][k]*seg[k]*σ[k,k]^3 for k in @groups(i)) for i in @comps)
    return val
end

function data(model::gcPCSAFTModel,V,T,z)
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

function a_hc(model::gcPCSAFTModel, V, T, z,_data=@f(data))
    ngroups = length(model.groups.flattenedgroups)

    _d,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    m = model.params.segment.values
    Σz = sum(z)
    c1 = 1/(1-ζ3)
    c2 = 3ζ2/(1-ζ3)^2
    c3 = 2ζ2^2/(1-ζ3)^3
    a_hs = bmcs_hs(ζ0,ζ1,ζ2,ζ3)
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

function a_disp(model::gcPCSAFTModel, V, T, z,_data=@f(data))
    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Σz = sum(z)
    I₁ = @f(I,1,_data)
    I₂ = @f(I,2,_data)
    C₁ = @f(C1,_data)
    m2ϵσ3₁,m2ϵσ3₂ = @f(m2ϵσ3)
    return -2*π*N_A*Σz/V*I₁*m2ϵσ3₁ - π*m̄*N_A*Σz/V*C₁*I₂*m2ϵσ3₂
end

function m2ϵσ3(model::gcPCSAFTModel, V, T, z)
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

function d(model::gcPCSAFTModel, V, T, z)
    ϵᵢᵢ = diagvalues(model.params.epsilon)
    σᵢᵢ = diagvalues(model.params.sigma) 
    return σᵢᵢ .* (1 .- 0.12 .* exp.(-3ϵᵢᵢ ./ T))
end

function ζ0123(model::gcPCSAFTModel, V, T, z,_d)
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

function Δ(model::gcPCSAFTModel, V, T, z, i, j, a, b,_data=@f(data))
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

function  Δ(model::gcPCSAFT, V, T, z,_data=@f(data))
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