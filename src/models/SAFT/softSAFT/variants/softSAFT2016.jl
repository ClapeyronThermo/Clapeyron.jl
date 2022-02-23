struct softSAFT2016Param <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type softSAFT2016Model <: softSAFTModel end

struct softSAFT2016{T} <: softSAFT2016Model
    components::Array{String,1}
    icomponents::UnitRange{Int}
    sites::SiteParam
    params::softSAFT2016Param
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
    lj::LJRefConsts
end

@registermodel softSAFT2016
function softSAFT2016(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())

    params,sites = getparams(components, ["SAFT/softSAFT"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["m"]
    k = params["k"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    packagedparams = softSAFT2016Param(segment, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1080/002689797170707"]
    icomponents = 1:length(components)
    return softSAFT2016(components,icomponents,sites,packagedparams,init_idealmodel,assoc_options,references, LJRefConsts())
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
