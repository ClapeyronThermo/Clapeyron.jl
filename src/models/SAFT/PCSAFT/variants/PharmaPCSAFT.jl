struct pharmaPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    k::PairParam{Float64}
    kT::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type pharmaPCSAFTModel <: PCSAFTModel end

struct pharmaPCSAFT{T <: IdealModel} <: pharmaPCSAFTModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    sites::SiteParam
    params::PCSAFTParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
    water::SpecialComp
end

@registermodel pharmaPCSAFT

export pharmaPCSAFT

function pharmaPCSAFT(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions(combining = :elliott))

    params,sites = getparams(components, ["SAFT/PCSAFT","properties/molarmass.csv"]; 
    userlocations=userlocations, 
    verbose=verbose,
    ignore_missing_singleparams = ["kT"])
    
    water = SpecialComp(components,["water08"])
    icomponents = 1:length(components)
    segment = params["m"]
    k0 = params["k"]
    n = length(components)
    k1 = get(params,"kT",PairParam("kT",components,zeros(n,n)))
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"])
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon,k0, k1, epsilon_assoc, bondvol)
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]

    model = PCSAFT(components,icomponents,sites,packagedparams,init_idealmodel,assoc_options,references,water)
    return model
end

#Δσh20(T) = σ[T] - σconstant
# σ = σconstant + Δσh20(T)
#https://doi.org/10.1016/j.cep.2007.02.034
Δσh20(T) = (10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10 
@inline water08_k(model::PCSAFTModel) = 0
@inline water08_k(model::PCSAFT) = model.water[]

function d(model::pharmaPCSAFTModel, V, T, z)
    ϵᵢᵢ = model.params.epsilon.diagvalues
    σᵢᵢ = model.params.sigma.diagvalues 
    _d = zeros(typeof(T),length(z))
    for i ∈ @comps
        _d[i] = σᵢᵢ[i]*(1- 0.12*exp(-3ϵᵢᵢ[i]/T))
    end
    k = water08_k(model)
    if !iszero(k)
        _d[k] += Δσh20(T)*(1- 0.12*exp(-3ϵᵢᵢ[k]/T))
    end
    return _d
end

function m2ϵσ3(model::pharmaPCSAFTModel, V, T, z)
    m = model.params.segment.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    k0 = model.params.k0.values
    k1 = model.params.k1.values
    m2ϵσ3₂ = zero(V+T+first(z))
    m2ϵσ3₁ = m2ϵσ3₂
    
    k = water08_k(model)

    @inbounds for i ∈ @comps
        zi = z[i]
        mi = m[i]
        ki = (k != i)
        constant = ki*zi*zi*mi*mi*(σ[i,i])^3
        exp1 = (ϵ[i,i]/T)
        exp2 = exp1*exp1
        m2ϵσ3₁ += constant*exp1
        m2ϵσ3₂ += constant*exp2
        for j ∈ 1:(i-1)
            kj = (j != k)
            constant = ki*kj*zi*z[j]*mi*m[j] * σ[i,j]^3 
            exp1 = ϵ[i,j]*(1 - k0[i,j] - k1[i,j]*T)/T
            exp2 = exp1*exp1
            m2ϵσ3₁ += 2*constant*exp1
            m2ϵσ3₂ += 2*constant*exp2
        end
    end

    if !iszero(k)
        Δσ = Δσh20(T)   
        zkmk = z[k]*m[k]
        @inbounds for j ∈ @comps
            σij = σ[k,j] + 0.5*((k==i) + (k==j))*Δσ
            constant = zkmk*z[j]*m[j]*σij^3
            exp1 = (ϵ[k,j]/T)
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
    Δσ = ifelse(iszero(k),zero(T),Δσh20(T))     
    σij = σ[i,j] + 0.5*((k==i) + (k==j))*Δσ
    res = gij*σij^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κijab
    return res
end

#optimized version, reduces memory allocations, but is specific to this PCSAFT
#thats why it is bound to the specific PCSAFT struct
#instead of the more general PCSAFTModel
function  Δ(model::pharmaPCSAFT, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    k = water08_k(model)
    k = model.water[]
    Δσ = ifelse(iszero(k),zero(T),Δσh20(T))     
    Δres = zero_assoc(κ,typeof(V+T+first(z)))
    for (idx,(i,j),(a,b)) in indices(Δres)
        gij = @f(g_hs,i,j,_data)
        σij = σ[i,j] + 0.5*((k==i) + (k==j))*Δσ
        Δres[idx] = gij*σij^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]
    end
    return Δres
end
