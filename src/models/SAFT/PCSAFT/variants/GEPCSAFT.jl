
abstract type GEPCSAFTModel <: SAFTModel end

const GEPCSAFTParam = PCSAFTParam

struct GEPCSAFT{I <: IdealModel,T,γ} <: GEPCSAFTModel
    components::Array{String,1}
    sites::SiteParam
    activity::γ
    params::PCSAFTParam{T}
    idealmodel::I
    assoc_options::AssocOptions
    references::Array{String,1}
end

"""
    GEPCSAFT <: SAFTModel

    GEPCSAFT(components;
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
- `activity`: Activity model

## Description

Perturbed-Chain SAFT (PC-SAFT), with Gᴱ mixing rule.

## References
1. Walker, P. J. (2022). Toward advanced, predictive mixing rules in SAFT equations of state. Industrial & Engineering Chemistry Research. [doi:10.1021/acs.iecr.2c03464](https://doi.org/10.1021/acs.iecr.2c03464)
"""
GEPCSAFT

export GEPCSAFT
function GEPCSAFT(components;
    idealmodel = BasicIdeal,
    activity = UNIFAC,
    userlocations = String[],
    ideal_userlocations = String[],
    activity_userlocations = String[],
    assoc_options = AssocOptions(),
    reference_state = nothing,
    verbose = false)

    params = getparams(components, ["SAFT/PCSAFT/PCSAFT_like.csv","SAFT/PCSAFT/PCSAFT_unlike.csv","SAFT/PCSAFT/PCSAFT_assoc.csv"]; userlocations = userlocations, verbose = verbose)
    sites = params["sites"]
    segment = params["segment"]
    k = get(params,"k",nothing)
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]

    packagedparams = GEPCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_activity = init_model(activity,components,activity_userlocations,verbose)
    references = String["10.1021/acs.iecr.2c03464"]
    model = GEPCSAFT(format_components(components),sites,init_activity,packagedparams,init_idealmodel,assoc_options,references)
    set_reference_state!(model,reference_state;verbose)
    return model
end

function _pcsaft(model::GEPCSAFT{I,T}) where {I,T}
    return PCSAFT{I,T}(model.components,model.sites,model.params,model.idealmodel,model.assoc_options,model.references)
end

function a_res(model::GEPCSAFTModel, V, T, z)    
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::GEPCSAFTModel,V,T,z)
    _d = @f(d)
    ζ0,ζ1,ζ2,ζ3 = @f(ζ0123,_d)
    m = model.params.segment.values
    m̄ = dot(z, m)/sum(z)
    return (_d,ζ0,ζ1,ζ2,ζ3,m̄)
end

function a_hc(model::GEPCSAFTModel, V, T, z , _data = @f(data))
    _,_,_,_,η,m̄ = _data
    g_hs = (1-η/2)/(1-η)^3
    a_hs = (4η-3η^2)/(1-η)^2
    return m̄*a_hs - (m̄-1)*log(g_hs)
end

function g_hs(model::GEPCSAFTModel, V, T, z,_data = @f(data))
    _,_,_,_,η,_ = _data
    return (1-η/2)/(1-η)^3
end

function a_hs(model::GEPCSAFTModel, V, T, z)
    _,_,_,_,η,_ = _data
    return (4η-3η^2)/(1-η)^2
end

function Δ(model::GEPCSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    g_hs_ = @f(g_hs,_data)
    return g_hs_*σij^3*(exp(ϵ_associjab/T)-1)*κijab
end

function a_disp(model::GEPCSAFTModel, V, T, z,_data=@f(data))
    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Σz = sum(z)
    m2ϵσ3₁,m2ϵσ3₂ = @f(m2ϵσ3,_data)
    πNAρ = π*N_A*Σz/V
    return -2*πNAρ*@f(I,1,_data)*m2ϵσ3₁ - m̄*πNAρ*@f(C1,_data)*@f(I,2,_data)*m2ϵσ3₂
end

d(model::GEPCSAFTModel, V, T, z) = ck_diameter(model, T, z)

function C1(model::GEPCSAFTModel, V, T, z, _data=@f(data))
    return C1(_pcsaft(model),V,T,z,_data)
end

function m2ϵσ3(model::GEPCSAFTModel, V, T, z,_data=@f(data))
    d,_,_,_,_,m̄ = _data
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    Σz = sum(z)
    _Ī = @f(Ī,1,m̄)
    #Iᵢ = @f(Ii,1,_data)
    md³ = zero(Base.promote_eltype(model,V,T,z))
    m²σ³ = zero(md³)
    A = zero(md³)
    @inbounds for i ∈ @comps
        mᵢ,dᵢ,σᵢ,ϵᵢ,zᵢ = m[i],d[i],σ[i,i],ϵ[i,i],z[i]
        σ³ᵢ,d³ᵢ = σᵢ*σᵢ*σᵢ,dᵢ*dᵢ*dᵢ
        m²σ³ += zᵢ*mᵢ*mᵢ*σ³ᵢ
        md³ += zᵢ*mᵢ*d³ᵢ
        A += zᵢ*mᵢ*σ³ᵢ/d³ᵢ*(ϵᵢ/T)*@f(Ī,1,mᵢ)
    end
    m²σ³,md³ = m²σ³/Σz,md³/Σz
    gₑ = excess_gibbs_free_energy(model.activity,V,T,z)/(R̄*T)
    m2ϵσ3₁ = md³/_Ī*(A - gₑ/12)/Σz
    m2ϵσ3₂ =  m2ϵσ3₁*m2ϵσ3₁/m²σ³
    return m2ϵσ3₁,m2ϵσ3₂
end

function I(model::GEPCSAFTModel, V, T, z, n, _data=@f(data))
    return I(_pcsaft(model),V,T,z,n,_data)
end

function Ī(model::GEPCSAFTModel, V, T, z, n, m)
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = 0.0*m
    @inbounds for i ∈ 1:7
        corr1,corr2,corr3 = corr[i]
        ki =corr1 + (m-1)/m*corr2 + (m-1)/m*(m-2)/m*corr3
        res += ki
    end
    return res
end
