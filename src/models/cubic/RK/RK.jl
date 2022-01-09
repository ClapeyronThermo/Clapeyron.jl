struct RKParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type RKModel <: ABCubicModel end

struct RK{T <: IdealModel,α,c,M} <: RKModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::M
    translation::c
    params::RKParam
    idealmodel::T
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel RK
export RK

function RK(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = RKAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(RK,init_mixing,Tc,pc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = RKParam(a,b,Tc,pc,Mw)
    references = String[]
    model = RK(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,1e-12,references)
    return model
end

function ab_consts(::Type{<:RKModel})
    Ωa =  1/(9*(2^(1/3)-1))
    Ωb = (2^(1/3)-1)/3
    return Ωa,Ωb
end

function cubic_abp(model::RKModel, V, T, z)
    n = sum(z)
    a,b,c = cubic_ab(model,V,T,z,n)
    v = V/n+c
    p =  R̄*T/(v-b) - a/((v+b)*v)
    return a,b,p
end

function cubic_poly(model::RKModel,p,T,z)
    n = sum(z)
    a,b,c = cubic_ab(model,p,T,z,n)
    RT⁻¹ = 1/(R̄*T)
    A = a*p* RT⁻¹* RT⁻¹
    B = b*p* RT⁻¹
    _1 = one(a)
    return (-A*B, -B*(B+_1) + A, -_1, _1),c
end


function a_res(model::RKModel, V, T, z)
    n=sum(z)
    ā,b̄,c̄ = cubic_ab(model,V,T,z,n)
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    RT⁻¹ = 1/(R̄*T)
    return -log(1+(c̄-b̄)*ρ) - ā*RT⁻¹*log(b̄*ρt+1)/b̄
    #return -log(V-n*b̄) - ā/(R̄*T*b̄*√(T/T̄c))*log(1+n*b̄/V)
end

cubic_zc(::RKModel) = 1/3

# include("variants/SRK.jl")
