abstract type ClausiusModel <: ABCCubicModel end

const ClausiusParam = ABCCubicParam

struct Clausius{T <: IdealModel,α,c,γ} <:ClausiusModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::γ
    translation::c
    params::ClausiusParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel Clausius
"""
    Clausius(components::Vector{String};
    idealmodel=BasicIdeal,
    alpha = NoAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=ClausiusTranslation,
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

Clausius Equation of state.
```
P = RT/(v-b) + a•α(T)/((v - Δ₀b)^2)

aᵢᵢ =27/64 * (RTcᵢ)²/Pcᵢ
bᵢᵢ = Vcᵢ - 1/4 * RTcᵢ/Pcᵢ
cᵢ = 3/8 * RTcᵢ/Pcᵢ - Vcᵢ

Δ₀ = ∑cᵢxᵢ/∑bᵢxᵢ

```

## References

1. Clausius, R. (1880). Ueber das Verhalten der Kohlensäure in Bezug auf Druck, Volumen und Temperatur. Annalen der Physik, 245(3), 337–357. [doi:10.1002/andp.18802450302](https://doi.org/10.1002/andp.18802450302)

"""
Clausius

export Clausius
function Clausius(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = ClausiusAlpha,
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
    a,b = ab_premixing(Clausius,init_mixing,Tc,pc,Vc,k)
    c = c_premixing(Clausius,init_mixing,Tc,pc,Vc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = ClausiusParam(a,b,c,Tc,pc,Vc,Mw)
    references = String["10.1002/andp.18802450302"]
    model = Clausius(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

function ab_premixing(model::Type{<:ClausiusModel},mixing::MixingRule,Tc,pc,vc,kij)
    _Tc = Tc.values
    _Vc = vc.values
    _pc = pc.values
    components = vc.components
    
    a = epsilon_LorentzBerthelot(SingleParam("a",components, @. 27*R̄^2*_Tc^2/_pc/64),kij)
    b = sigma_LorentzBerthelot(SingleParam("b",components, @. _Vc-1/4*R̄*_Tc/_pc))
    #a = epsilon_LorentzBerthelot(SingleParam("a",components, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    #b = sigma_LorentzBerthelot(SingleParam("b",components, @. Ωb*R̄*_Tc/_pc))
    return a,b
end

function c_premixing(model::Type{<:ClausiusModel},mixing::MixingRule,Tc,pc,vc,kij)
    _Tc = Tc.values
    _Vc = vc.values
    _pc = pc.values
    components = vc.components
 
    c = sigma_LorentzBerthelot(SingleParam("c",components, @. 3/8*R̄*_Tc/_pc-_Vc))
    #a = epsilon_LorentzBerthelot(SingleParam("a",components, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    #b = sigma_LorentzBerthelot(SingleParam("b",components, @. Ωb*R̄*_Tc/_pc))
    return c
end

function cubic_Δ(model::ClausiusModel,z) 
    b = model.params.b.diagvalues
    c = model.params.c.diagvalues
    z⁻¹ = sum(z)^-1
    b̄ = dot(b,z)*z⁻¹
    c̄ = dot(c,z)*z⁻¹
    return (-c̄/b̄,-c̄/b̄)
end

crit_pure(model::ClausiusModel) = crit_pure_tp(model)
#=
 (-B2-2(B2+B)+A)
 (-B2-2B2-2B+A)
 (-3B2-2B+A)
=#