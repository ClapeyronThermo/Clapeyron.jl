


struct PRParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type PRModel <: ABCubicModel end

struct PR{T <: IdealModel,α,γ} <: RKModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    activity::γ
    params::PRParam
    idealmodel::T
    absolutetolerance::Float64
    references::Array{String,1}
end

has_sites(::Type{<:PRModel}) = false
has_groups(::Type{<:PRModel}) = false
built_by_macro(::Type{<:PRModel}) = false

function Base.show(io::IO, mime::MIME"text/plain", model::PR)
    return eosshow(io, mime, model)
end

function Base.show(io::IO, model::PR)
    return eosshow(io, model)
end

Base.length(model::PR) = Base.length(model.icomponents)

molecular_weight(model::PR,z=SA[1.0]) = group_molecular_weight(model.groups,mw(model),z)

export PR
function PR(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = SoaveAlpha,
    activity = nothing,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    activity_userlocations = String[],
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    _pc = params["pc"]
    pc = _pc.values
    Mw = params["Mw"]
    _Tc = params["Tc"]
    Tc = _Tc.values
    Ωa, Ωb = ab_consts(PR)
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. Ωa*R̄^2*Tc^2/pc, k))
    #check if this is correct in the general case.
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. Ωb*R̄*Tc/pc))
    
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_activity = init_model(activity,components,activity_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = PRParam(a, b, params["Tc"],_pc,Mw)
    references = String[]
    model = PR(components,icomponents,init_alpha,init_activity,packagedparams,init_idealmodel,1e-12,references)
    return model
end


function ab_consts(::Type{<:PRModel})
    return 0.457235,0.077796
end

function cubic_ab(model::PR{<:Any,SoaveAlpha},T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    α = model.alpha
    ω = α.params.acentricfactor.values
    Tc = model.params.Tc.values
    invn = one(n)/n
    αx = @. (1+(0.37464+1.54226*ω-0.26992*ω^2)*(1-√(T/Tc))) * z * invn 
    āᾱ = dot(αx,Symmetric(a),αx)
    b̄ = dot(z,Symmetric(b),z) * invn*invn
    return āᾱ ,b̄



function cubic_abp(model::PRModel, V, T, z) 
    n = sum(z)
    v = V/n
    āᾱ ,b̄ = cubic_ab(model,T,z,n)
    _1 = one(b̄)
    denom = evalpoly(v,(-b̄*b̄,2*b̄,_1))
    p = R̄*T/(v-b̄) - āᾱ /denom
    return āᾱ, b̄,p
end

function cubic_poly(model::PRModel,p,T,z)
    a,b = cubic_ab(model,T,z)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    k₀ = B*(B*(B+1.0)-A)
    k₁ = -B*(3*B+2.0) + A
    k₂ = B-1.0
    k₃ = 1.0
    return [k₀,k₁,k₂,k₃]
end

#=
 (-B2-2(B2+B)+A)
 (-B2-2B2-2B+A)
 (-3B2-2B+A)
=#
function a_res(model::PRModel, V, T, z)
    n = sum(z)
    a,b = cubic_ab(model,T,z,n)
    Δ1 = 1+√2
    Δ2 = 1-√2
    ΔPRΔ = 2*√2
    RT⁻¹ = 1/(R̄*T)
    ρ = n/V
    return -log(1-b*ρ) - a*RT⁻¹*log((Δ1*b*ρ+1)/(Δ2*b*ρ+1))/(ΔPRΔ*b)

    #return -log(V-n*b̄) + āᾱ/(R̄*T*b̄*2^(3/2)) * log((2*V-2^(3/2)*b̄*n+2*b̄*n)/(2*V+2^(3/2)*b̄*n+2*b̄*n))
end

cubic_zc(::PRModel) = 0.3074
