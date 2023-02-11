struct softSAFT2016Param <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type softSAFT2016Model <: softSAFTModel end

struct softSAFT2016{T} <: softSAFT2016Model
    components::Array{String,1}
    sites::SiteParam
    params::softSAFT2016Param
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
    lj::LJRefConsts
end

@registermodel softSAFT2016

"""
    softSAFT2016Model <: softSAFTModel

    softSAFT2016(components; 
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

Soft SAFT, with Lennard-Jones function from Thol et al. (2016)

## References
1. FELIPE J. BLAS and LOURDES F. VEGA. (1997). Thermodynamic behaviour of homonuclear and heteronuclear Lennard-Jones chains with association sites from simulation and theory. Molecular physics, 92(1), 135–150. [doi:10.1080/002689797170707](https://doi.org/10.1080/002689797170707)
2. Thol, M., Rutkai, G., Köster, A., Lustig, R., Span, R., & Vrabec, J. (2016). Equation of state for the Lennard-Jones fluid. Journal of physical and chemical reference data, 45(2), 023101. [doi:10.1063/1.4945000](https://doi.org/10.1063/1.4945000)
"""
softSAFT2016

function softSAFT2016(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions(), kwargs...)

    params,sites = getparams(components, ["SAFT/softSAFT","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["segment"]
    k = get(params,"k",nothing)
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options) #combining rules for association
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    packagedparams = softSAFT2016Param(params["Mw"],segment, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1080/002689797170707","10.1063/1.4945000"]
    return softSAFT2016(components,sites,packagedparams,init_idealmodel,assoc_options,references, LJRefConsts())
end

function a_LJ(model::softSAFT2016Model,V,T,z,_data = @f(data))
    σ3,ϵ,m̄,_ = _data
    V0 = m̄*N_A*σ3
    T0 = ϵ
    τ = 1.32/(T/T0)
    Vx = V/sum(z)
    δ = (V0/Vx)/0.31 
    ai = zero(δ+τ)
    n = model.lj.n
    t = model.lj.t
    d = model.lj.d
    c = model.lj.c
    β = model.lj.beta
    γ = model.lj.gamma
    η = model.lj.eta
    ε = model.lj.epsilon   

    @inbounds begin
        for k ∈ 1:6
            ai += n[k]*(δ^d[k])*(τ^t[k])
        end
        for (k,k_) ∈ zip(7:12,1:6)
            ai += n[k]*(δ^d[k])*(τ^t[k])*exp(-δ^c[k_])
        end

        for (k,k_) ∈ zip(13:23,1:11)
            ai += n[k]*(δ^(d[k]))*(τ^(t[k]))*
            exp(-η[k_]*(δ - ε[k_])^2 - β[k_]*(τ -γ[k_])^2)
        end
    end
    return m̄*ai
end

export softSAFT2016
