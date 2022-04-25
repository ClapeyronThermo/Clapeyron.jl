
abstract type PRModel <: ABCubicModel end

const PRParam = ABCubicParam

struct PR{T <: IdealModel,α,c,γ} <:PRModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::γ
    translation::c
    params::PRParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel PR

"""
    PR(components::Vector{String}; idealmodel=BasicIdeal,
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
Peng-Robinson Equation of state.
```
P = RT/(V-Nb) + a•α(T)/(V-Nb₁)(V-Nb₂)
b₁ = (1 + √2)b
b₂ = (1 - √2)b
```

## References
1. Peng, D.Y., & Robinson, D.B. (1976). A New Two-Constant Equation of State. Industrial & Engineering Chemistry Fundamentals, 15, 59-64. doi:10.1021/I160057A011
"""
PR


export PR
function PR(components::Vector{String}; idealmodel=BasicIdeal,
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
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(PR,init_mixing,Tc,pc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = PRParam(a,b,Tc,pc,Mw)
    references = String["10.1021/I160057A011"]
    model = PR(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

function ab_consts(::Type{<:PRModel})
    return 0.457235,0.077796
end

function cubic_Δ(model::PRModel,z) 
    sqrt2 = sqrt(2)
    return (1+sqrt2,1-sqrt2)
end

#=
 (-B2-2(B2+B)+A)
 (-B2-2B2-2B+A)
 (-3B2-2B+A)
=#