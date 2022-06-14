const RKParam = ABCubicParam
abstract type RKModel <: ABCubicModel end

struct RK{T <: IdealModel,α,c,M} <: RKModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::M
    translation::c
    params::RKParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel RK
export RK

"""
    RK(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = PRAlpha,
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

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)

## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description
Redlich-Kwong Equation of state.
```
P = RT/(V-Nb) + a•α(T)/(V(V+Nb))
```

## References
1. Redlich, O., & Kwong, J. N. S. (1949). On the thermodynamics of solutions; an equation of state; fugacities of gaseous solutions. Chemical Reviews, 44(1), 233–244. doi:10.1021/cr60137a013
"""
RK

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
    references = String["10.1021/cr60137a013"]
    model = RK(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
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
    _1 = one(A)
    return (-A*B, -B*(B+_1) + A, -_1, _1),c
end


function a_res(model::RKModel, V, T, z,_data = data(model,V,T,z))
    n,ā,b̄,c̄ = _data
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    RT⁻¹ = 1/(R̄*T)
    return -log(1+(c̄-b̄)*ρ) - ā*RT⁻¹*log(b̄*ρt+1)/b̄
    #return -log(V-n*b̄) - ā/(R̄*T*b̄*√(T/T̄c))*log(1+n*b̄/V)
end

cubic_zc(::RKModel) = 1/3

# include("variants/SRK.jl")
