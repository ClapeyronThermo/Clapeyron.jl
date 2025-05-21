"""
    sCPAModel <: CPAModel

    function sCPA(components;
        idealmodel = BasicIdeal,
        radial_dist::Symbol = :KG,
        cubicmodel = RK,
        alpha = sCPAAlpha,
        mixing = vdW1fRule,
        activity = nothing,
        translation = NoTranslation,
        userlocations = String[],
        ideal_userlocations = String[],
        alpha_userlocations = String[],
        activity_userlocations = String[],
        mixing_userlocations = String[],
        translation_userlocations = String[],
        reference_state = nothing,
        verbose = false,
        assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `a`: Single Parameter (`Float64`) - Atraction parameter `[m^6*Pa/mol]`
- `b`: Single Parameter (`Float64`) - Covolume `[m^3/mol]`
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`) - Mixed Atraction Parameter `[m^6*Pa/mol]`
- `b`: Pair Parameter (`Float64`) - Mixed Covolume `[m^3/mol]`
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[J]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Input models
- `idealmodel`: Ideal Model
- `cubicmodel`: Cubic Model

## Description

Simplified Cubic Plus Association (s-CPA) EoS. Consists in the addition of a cubic part and an association part:
```
a_res(model::CPA) = a_res(model::Cubic) + a_assoc(model)
```

The `radial_dist` argument can be used to choose between a Carnahan-Starling form (`CS`, default) or the Kontogeorgis (`KG`) term, more widely known as s-CPA. using `sCPA(components, radial_dist =: CS)` is equivalent to using the original CPA

## References
1. Kontogeorgis, G. M., Michelsen, M. L., Folas, G. K., Derawi, S., von Solms, N., & Stenby, E. H. (2006). Ten years with the CPA (cubic-plus-association) equation of state. Part 1. Pure compounds and self-associating systems. Industrial & Engineering Chemistry Research, 45(14), 4855–4868. [doi:10.1021/ie051305v](https://doi.org/10.1021/ie051305v)
"""
function sCPA(components;
            idealmodel = BasicIdeal,
            radial_dist::Symbol = :KG,
            cubicmodel = RK,
            alpha = sCPAAlpha,
            mixing = vdW1fRule,
            activity = nothing,
            translation = NoTranslation,
            userlocations = String[],
            ideal_userlocations = String[],
            alpha_userlocations = String[],
            activity_userlocations = String[],
            mixing_userlocations = String[],
            translation_userlocations = String[],
            reference_state = nothing,
            verbose = false,
            assoc_options = AssocOptions())

    return CPA(components;
        idealmodel = idealmodel,
        radial_dist = radial_dist,
        cubicmodel = cubicmodel,
        alpha = alpha,
        mixing = mixing,
        activity = activity,
        translation = translation,
        userlocations = userlocations,
        ideal_userlocations = ideal_userlocations,
        alpha_userlocations = alpha_userlocations,
        activity_userlocations = activity_userlocations,
        mixing_userlocations = mixing_userlocations,
        translation_userlocations = translation_userlocations,
        reference_state = reference_state,
        verbose = verbose,
        assoc_options = assoc_options)
end

default_locations(::typeof(sCPA)) = ["SAFT/CPA/sCPA/", "properties/molarmass.csv","properties/critical.csv"]

export sCPA
