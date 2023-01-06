abstract type KUModel <: ABCubicModel end

struct KUParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    omega_a::SingleParam{Float64}
    omega_b::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct KU{T <: IdealModel,α,c,γ} <:KUModel
    components::Array{String,1}
    alpha::α
    mixing::γ
    translation::c
    params::KUParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel KU

"""
    KU(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = KUAlpha,
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
- `vc`: Single Parameter (`Float64`) - Critical Volume `[m^3]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`) (optional)
- `l`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m^3]`
- `omega_a`: Single Parameter (`Float64`) - Critical Constant for a - No units
- `omega_b`: Single Parameter (`Float64`) - Critical Constant for b - No units
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
Kumar-Upadhyay Cubic Equation of state. `Ωa` and `Ωb` are component-dependent
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
1. Kumar, A., & Upadhyay, R. (2021). A new two-parameters cubic equation of state with benefits of three-parameters. Chemical Engineering Science, 229(116045), 116045. [doi:10.1016/j.ces.2020.116045](https://doi.org/10.1016/j.ces.2020.116045)
"""
KU
export KU

#another alternative would be to store the Ωa, Ωb in the mixing struct.

function KU(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = KUAlpha,
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
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    Vc = params["Vc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    
    n = length(components)
    a = PairParam("a",components,zeros(n))
    b = PairParam("b",components,zeros(n))
    omega_a = SingleParam("Ωa",components,zeros(n))
    omega_b = SingleParam("Ωb",components,zeros(n))
    
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    packagedparams = KUParam(a,b,omega_a,omega_b,Tc,pc,Vc,Mw)
    references = String["10.1016/j.ces.2020.116045"]
    model = KU(components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    return model
end

function ab_premixing(model::KUModel,mixing::MixingRule,k,l)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _vc = model.params.Vc
    a = model.params.a
    b = model.params.b
    Zc = _pc .* _vc ./ (R̄ .* _Tc)
    χ  = @. cbrt(sqrt(1458*Zc^3 - 1701*Zc^2 + 540*Zc -20)/(32*sqrt(3)*Zc^2) 
    - (729*Zc^3 - 216*Zc + 8)/(1728*Zc^3))
    α  = @. (χ + (81*Zc^2 - 72*Zc + 4)/(144*χ*Zc^2) + (3*Zc - 2)/(12*Zc))
    Ωa = @. Zc*(1 + 1.6*α - 0.8*α^2)^2/(1 - α)^2/(2 + 1.6*α)
    Ωb = @. Zc*α
    model.params.omega_a .= Ωa
    model.params.omega_b .= Ωb
    diagvalues(a) .= @. Ωa*R̄^2*_Tc^2/_pc
    diagvalues(b) .= @. Ωb*R̄*_Tc/_pc
    a = epsilon_LorentzBerthelot!(a,k)
    b = sigma_LorentzBerthelot!(b,l)
    return a,b
end

ab_consts(model::KUModel) = model.params.omega_a.values,model.params.omega_b.values

#only used in premixing
cubic_Δ(model::KUModel,z) = (0.4,-2.0)

function T_scale(model::KUModel,z=SA[1.0])
    Tc,_ = vdw_tv_mix(model.params.Tc.values,model.params.Vc.values,z)
    return Tc
end

function p_scale(model::KUModel,z=SA[1.0])
    return dot(model.params.Pc.values,z)/sum(z)
end

kumar_zc(model::KUModel) = only(model.params.Pc.values)*only(model.params.Vc.values)/(R̄*only(model.params.Tc.values))

function x0_crit_pure(model::KUModel)
    lb_v = lb_volume(model)
    vc = model.params.Vc.values[1]
    (1.1, log10(vc))
end