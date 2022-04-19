abstract type KumarCubicModel <: ABCubicModel end

struct KumarCubicParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct KumarCubic{T <: IdealModel,α,c,γ} <:KumarCubicModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::γ
    translation::c
    params::KumarCubicParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel KumarCubic

"""
    KumarCubic(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = KCAlpha,
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
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m^3]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m^3]`
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
Ashutosh Kumar and Rajeev Upadhyay "Two-Parameter Cubic" Equation of state.
```
P = RT/(v-b) + a•κ(T)/((v²-1.6bv - 0.8b²)
a = Ωa(R²Tcᵢ²/Pcᵢ)
b = Ωb(R²Tcᵢ²/Pcᵢ)
Ωa = Zc[(1 + 1.6α - 0.8α²)²/((1 - α²)(2 + 1.6α))]
χ = ∛[√(1458Zc³ - 1701Zc² + 540Zc -20)/32√3Zc² - (729Zc³ - 216Zc + 8)/1728Zc³]
α = [χ + (81Zc² - 72Zc + 4)/144Zc²χ + (3Zc - 2)/12Zc]
Ωa = Zc[(1 + 1.6α - 0.8α²)²/((1 - α²)(2 + 1.6α))]
Ωb = αZc
```

## References
1. Kumar, A., & Upadhyay, R. (2021). A new two-parameters cubic equation of state with benefits of three-parameters. Chemical Engineering Science, 229(116045), 116045. doi:10.1016/j.ces.2020.116045
"""
KumarCubic

#another alternative would be to store the Ωa, Ωb in the mixing struct.
function ab_premixing(::Type{<:KumarCubic},mixing,Tc,pc,kij,Vc)
    Ωa, Ωb = ab_consts(T)
    _Tc = Tc.values
    _pc = pc.values
    _vc = Vc.values
    Zc = _pc .* _vc ./ (R̄ .* _tc)
    χ  = @. cbrt(sqrt(1458*Zc^3 - 1701*Zc^2 + 540*Z -20)/(32*sqrt(3)*Zc^2) - (729*Zc^3 - 216*Zc + 8)/(1728*Zc^3))
    α  = @. (χ + (81*Zc^2 - 72*Zc + 4)/(144*χ*Zc^2) + (3*Zc - 2)/(12*Zc))
    Ωa = @. Zc*((1 + 1.6*α - 0.8*α^2)^2/((1 - α^2)(2 + 1.6*α)))
    Ωb = @. Zc*α
    a = epsilon_LorentzBerthelot(SingleParam(pc, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    b = sigma_LorentzBerthelot(SingleParam(pc, @. Ωb*R̄*_Tc/_pc))
    return a,b
end

export KumarCubic
function KumarCubic(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = KCAlpha,
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
    Vc = params["Vc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(KumarCubic,init_mixing,Tc,pc,k,Vc)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = KumarCubicParam(a,b,Tc,pc,Vc,Mw)
    references = String["10.1016/j.ces.2020.116045"]
    model = KumarCubic(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

#only used in premixing
function ab_consts(::Type{<:KumarCubicModel})
    return 1.0,1.0
end

cubic_Δ(model::KumarCubicModel) = (-0.4,2.0)

#only used for single component properties, we need a better abstraction here
cubic_zc(model::KumarCubicModel) = only(model.Pc.values)*only(model.Vc.values)/(R̄*only(model.Tc.values))
