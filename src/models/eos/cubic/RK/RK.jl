struct RKParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type RKModel <: ABCubicModel end

struct RK{T <: IdealModel,α,γ} <: RKModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::γ
    params::RKParam
    idealmodel::T
    absolutetolerance::Float64
    references::Array{String,1}
end

has_sites(::Type{<:RKModel}) = false
has_groups(::Type{<:RKModel}) = false
built_by_macro(::Type{<:RKModel}) = false

function Base.show(io::IO, mime::MIME"text/plain", model::RK)
    return eosshow(io, mime, model)
end

function Base.show(io::IO, model::RK)
    return eosshow(io, model)
end

Base.length(model::RK) = Base.length(model.icomponents)

molecular_weight(model::RK,z=SA[1.0]) = group_molecular_weight(model.groups,mw(model),z)

export RK
function RK(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = RKAlpha,
    mixing = vdW1fRule,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    _pc = params["pc"]
    pc = _pc.values
    Mw = params["Mw"]
    _Tc = params["Tc"]
    Tc = _Tc.values
    #T̄c = sum(sqrt.(Tc*Tc')) #is this term correctly calculated? sqrt(Tc*Tc') is a matrix sqrt
    Ωa, Ωb = ab_consts(RK)
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. Ωa*R̄^2*Tc^2/pc), k) 
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. Ωb*R̄*Tc/pc))
    
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_mixing = init_model(mixing,components,mixing_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = RKParam(a, b, params["Tc"],_pc,Mw)
    references = String[]
    model = RK(components,icomponents,init_alpha,init_mixing,packagedparams,init_idealmodel,1e-12,references)
    return model
end

function ab_consts(::Type{<:RKModel})
    Ωa =  1/(9*(2^(1/3)-1))
    Ωb = (2^(1/3)-1)/3
    return Ωa,Ωb
end
function cubic_ab(model::RK{<:Any,<:Any,<:Any},V,T,z=SA[1.0],n=sum(z))
    invn2 = (one(n)/n)^2
    a = model.params.a.values
    b = model.params.b.values
    α = @f(α_function,model.alpha)
    ā,b̄ = @f(mixing_rule,model.mixing,α,a,b)
    return ā ,b̄
end

function cubic_abp(model::RKModel, V, T, z)
    n = sum(z)
    a,b = cubic_ab(model,V,T,z,n)
    v = V/n
    p =  R̄*T/(v-b) - a/((v+b)*v)
    return a,b,p
end

function cubic_poly(model::RKModel,p,T,z)
    n = sum(z)
    a,b = cubic_ab(model,p,T,z,n)
    RT⁻¹ = 1/(R̄*T)
    A = a*p* RT⁻¹* RT⁻¹
    B = b*p* RT⁻¹
    _1 = one(a)
    return [-A*B, -B*(B+_1) + A, -_1, _1]
end


function a_res(model::RKModel, V, T, z)
    n = sum(z)
    ā,b̄ = cubic_ab(model,T,z,n)
    ρ = n/V
    RT⁻¹ = 1/(R̄*T)
    return -log(1-b̄*ρ) - ā*RT⁻¹*log(b̄*ρ+1)/b̄
    #return -log(V-n*b̄) - ā/(R̄*T*b̄*√(T/T̄c))*log(1+n*b̄/V)
end

cubic_zc(::RKModel) = 1/3

# include("variants/SRK.jl")
