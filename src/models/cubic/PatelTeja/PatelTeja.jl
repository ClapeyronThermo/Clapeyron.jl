abstract type PatelTejaModel <: ABCCubicModel end

const PatelTejaParam = ABCCubicParam

struct PatelTeja{T <: IdealModel,α,c,γ} <:PatelTejaModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::γ
    translation::c
    params::PatelTejaParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel PatelTeja
"""
    PatelTeja(components::Vector{String};
    idealmodel=BasicIdeal,
    alpha = NoAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=PatelTejaTranslation,
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

## Description

PatelTeja Equation of state. it uses [`vdW`](@ref) the following models:
- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`PatelTejaAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)
```
P = RT/(V-Nb) + a•α(T)/(V-c)²
```

## References

1. PatelTeja, D. (1899). Sur une méthode purement physique pour la détermination des poids moléculaires des gaz et des poids atomiques de leurs éléments. Journal de Physique Théorique et Appliquée, 8(1), 263–274. doi:10.1051/jphystap:018990080026300

"""
PatelTeja

export PatelTeja
function PatelTeja(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = PatelTejaAlpha,
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
    Vc = params["vc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(PatelTeja,init_mixing,Tc,pc,Vc,k)
    c = c_premixing(PatelTeja,init_mixing,Tc,pc,Vc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = PatelTejaParam(a,b,c,Tc,pc,Vc,Mw)
    references = String["10.1016/0009-2509(82)80099-7"]
    model = PatelTeja(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

function ab_premixing(model::Type{<:PatelTejaModel},mixing::MixingRule,Tc,pc,vc,kij)
    _Tc = Tc.values
    _Vc = vc.values
    _pc = pc.values
    components = vc.components
    _Zc = _pc.*_Vc./(R̄*_Tc)
             
    _poly = [(-_Zc[i]^3,3*_Zc[i]^2,2-3*_Zc[i],1.) for i ∈ 1:length(components)]
    
    sols = Solvers.roots3.(_poly)

    Ωb = [minimum(real.(sols[i][isreal.(sols[1]).*real.(sols[1]).>0])) for i ∈ 1:length(components)]
    Ωa = @. 3*_Zc^2+3*(1-2*_Zc)*Ωb+Ωb^2+1-3*_Zc
    
    a = epsilon_LorentzBerthelot(SingleParam("a",components, @. R̄^2*_Tc^2/_pc*Ωa),kij)
    b = sigma_LorentzBerthelot(SingleParam("b",components, @. R̄*_Tc/_pc*Ωb))
    return a,b
end

function c_premixing(model::Type{<:PatelTejaModel},mixing::MixingRule,Tc,pc,vc,kij)
    _Tc = Tc.values
    _Vc = vc.values
    _pc = pc.values
    components = vc.components
    _Zc = _pc.*_Vc./(R̄*_Tc)

    Ωc = @. 1-3*_Zc

    c = sigma_LorentzBerthelot(SingleParam("c",components, @. Ωc*R̄*_Tc/_pc))
    return c
end

function cubic_Δ(model::PatelTejaModel,z) 
    b = model.params.b.values
    c = model.params.c.values
    z⁻¹ = sum(z)^-1
    b̄ = ∑(b.*z)*z⁻¹
    c̄ = ∑(c.*z)*z⁻¹

    cb⁻¹ = c̄/b̄

    det12 = 1+6*cb⁻¹+cb⁻¹^2
    tr12 = 1+cb⁻¹
    return (-1/2*(tr12+sqrt(det12)),-1/2*(tr12-sqrt(det12)))
end

crit_pure(model::PatelTejaModel) = crit_pure_tp(model)
#=
 (-B2-2(B2+B)+A)
 (-B2-2B2-2B+A)
 (-3B2-2B+A)
=#