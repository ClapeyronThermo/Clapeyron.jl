
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

Clausius Equation of state. it uses [`vdW`](@ref) the following models:
- Translation Model: [`ClausiusTranslation`](@ref)
- Alpha Model: [`ClausiusAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)
```
P = RT/(V-Nb) + a•α(T)/(V-c)²
```

## References

1. Clausius, D. (1899). Sur une méthode purement physique pour la détermination des poids moléculaires des gaz et des poids atomiques de leurs éléments. Journal de Physique Théorique et Appliquée, 8(1), 263–274. doi:10.1051/jphystap:018990080026300

"""

function Clausius(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = ClausiusAlpha,
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
    
    return vdW(components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing=mixing,
    activity = activity,
    translation=translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    translation_userlocations = translation_userlocations,
    verbose = verbose)
end

export Clausius


