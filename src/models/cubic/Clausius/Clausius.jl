abstract type ClausiusModel <: ABCCubicModel end

struct Clausius{T <: IdealModel,α,c,M} <: ClausiusModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::M
    translation::c
    params::ABCCubicParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel Clausius
export Clausius

"""
    Clausius(components::Vector{String};
    idealmodel=BasicIdeal,
    alpha = NoAlpha,
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
- `vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
- `k`: Pair Parameter (`Float64`)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)
- `c`: Pair Parameter (`Float64`)

## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description

Clausius Equation of state.

```
P = RT/(V-Nb) + a•α(T)/V²
```

## References

1. Clausius, D. (1899). Sur une méthode purement physique pour la détermination des poids moléculaires des gaz et des poids atomiques de leurs éléments. Journal de Physique Théorique et Appliquée, 8(1), 263–274. doi:10.1051/jphystap:018990080026300

"""
Clausius

function Clausius(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = ClausiusAlpha,
    mixing = ClausiusRule,
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
    Vc = params["vc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(Clausius,init_mixing,Tc,pc,Vc,k)
    c = c_premixing(Clausius,Tc,Pc,Vc)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = ABCCubicParam(a,b,c,Tc,pc,Vc,Mw)
    references = String["10.1051/jphystap:018990080026300"]
    model = Clausius(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

function ab_premixing(model::Type{<:ClausiusModel},mixing::MixingRule,Tc,pc,vc,kij)
    _Tc = Tc.values
    _pc = pc.values
    _Vc = Vc.values
    components = Vc.components
    a = epsilon_LorentzBerthelot(SingleParam("a",components, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    b = epsilon_LorentzBerthelot(SingleParam("b",components, @. _Vc - 0.25*R̄*_Tc/_pc),kij)
    return a,b
end

function c_premixing(model::Type{<:ClausiusModel},Tc,Pc,vc)
    c = @* (3/8)*R̄*Tc/Pc - vc
end

function pure_cubic_zc(model::ClausiusModel)
    return only(model.params.Pc.values)*only(model.params.Vc.values)/(R̄*only(model.params.Tc.values))
end

function cubic_Δ(model,z)
    c = dot(model.params.c.values,z)/sum(z)
    return (c,c)
end