abstract type KUModel <: ABCCubicModel end

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

"""
    KU(components;
    idealmodel = BasicIdeal,
    alpha = KUAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional)
- `l`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³]`
- `omega_a`: Single Parameter (`Float64`) - Critical Constant for a - No units
- `omega_b`: Single Parameter (`Float64`) - Critical Constant for b - No units
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
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

## Model Construction Examples
```julia
# Using the default database
model = KU("water") #single input
model = KU(["water","ethanol"]) #multiple components
model = KU(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = KU(["water","ethanol"],alpha = TwuAlpha) #modifying alpha function
model = KU(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = KU(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = KU(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = PR78Alpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model = KU(["ethane","butane"],alpha = my_alpha)

# User-provided parameters, passing files or folders
model = KU(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = KU(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Vc = [4.25e-5, 6.43e-5],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.], #k,l can be ommited in single-component models.
                        l = [0. 0.01; 0.01 0.])
                    )
```

## References
1. Kumar, A., & Upadhyay, R. (2021). A new two-parameters cubic equation of state with benefits of three-parameters. Chemical Engineering Science, 229(116045), 116045. [doi:10.1016/j.ces.2020.116045](https://doi.org/10.1016/j.ces.2020.116045)
"""
KU
export KU

#another alternative would be to store the Ωa, Ωb in the mixing struct.

function KU(components;
    idealmodel = BasicIdeal,
    alpha = KUAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    formatted_components = format_components(components)
    params = getparams(formatted_components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations = userlocations, verbose = verbose)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    Vc = params["Vc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    n = length(formatted_components)
    a = PairParam("a",formatted_components,zeros(n,n),ones(Bool,n,n))
    b = PairParam("b",formatted_components,zeros(n,n),ones(Bool,n,n))
    omega_a = SingleParam("Ωa",formatted_components,zeros(n))
    omega_b = SingleParam("Ωb",formatted_components,zeros(n))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_alphamodel(alpha,components,params,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    packagedparams = KUParam(a,b,omega_a,omega_b,Tc,pc,Vc,Mw)
    references = String["10.1016/j.ces.2020.116045"]
    model = KU(formatted_components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
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
    ab_diagvalues!(a,b,Ωa,Ωb,_Tc,_pc,Rgas(model))
    a = epsilon_LorentzBerthelot!(a,k)
    b = sigma_LorentzBerthelot!(b,l)
    return a,b
end

function recombine_mixing!(model::KUModel,mixing_model,k = nothing,l = nothing)
    recombine!(mixing_model)
    a,b = ab_premixing(model,mixing_model,k,l)
    #we set this again just in case
    model.params.a .= a
    model.params.b .= b
    return mixing_model
end

ab_consts(model::KUModel) = model.params.omega_a.values,model.params.omega_b.values

#only used in premixing
cubic_Δ(model::KUModel,z) = (0.4,-2.0)

function T_scale(model::KUModel,z)
    Tc,_ = vdw_tv_mix(model.params.Tc.values,model.params.Vc.values,z)
    return Tc
end

function p_scale(model::KUModel,z)
    return dot(model.params.Pc.values,z)/sum(z)
end

function x0_crit_pure(model::KUModel,z)
    (1.1, log10(dot(model.params.Vc.values,z)/sum(z)))
end