abstract type ClausiusModel <: ABCCubicModel end

const ClausiusParam = ABCCubicParam

struct Clausius{T <: IdealModel,α,c,γ} <:ClausiusModel
    components::Array{String,1}
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
- `Vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`) (optional)
- `l`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
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
    k  = get(params,"k",nothing)
    l  = get(params,"l",nothing)
    pc = params["Pc"]
    Vc = params["Vc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",components,zeros(length(components)))
    b = PairParam("b",components,zeros(length(components)))
    c = PairParam("c",components,zeros(length(components)))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    packagedparams = ClausiusParam(a,b,c,Tc,pc,Vc,Mw)
    references = String["10.1002/andp.18802450302"]
    model = Clausius(components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    return model
end

function ab_premixing(model::ClausiusModel,mixing::MixingRule,k=nothing,l=nothing)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    a = model.params.a
    b = model.params.b

    diagvalues(a) .= @. 27*R̄^2*_Tc^2/_pc/64
    diagvalues(b) .= @. _Vc-1/4*R̄*_Tc/_pc
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    #a = epsilon_LorentzBerthelot(SingleParam("a",components, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    #b = sigma_LorentzBerthelot(SingleParam("b",components, @. Ωb*R̄*_Tc/_pc))
    return a,b
end

function c_premixing(model::ClausiusModel)
    _Tc = model.params.Tc
    _Vc = model.params.Vc
    _pc = model.params.Pc
    c = model.params.c
    diagvalues(c) .=  @. 3/8*R̄*_Tc/_pc-_Vc
    sigma_LorentzBerthelot!(c)
    #a = epsilon_LorentzBerthelot(SingleParam("a",components, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    #b = sigma_LorentzBerthelot(SingleParam("b",components, @. Ωb*R̄*_Tc/_pc))
    return c
end

function cubic_Δ(model::ClausiusModel,z) 
    b = diagvalues(model.params.b)
    c = diagvalues(model.params.c)
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