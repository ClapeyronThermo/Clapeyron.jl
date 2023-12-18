"""
    sCPAModel <: CPAModel

    function sCPA(components;
        idealmodel=BasicIdeal,
        cubicmodel=RK,
        alpha=sCPAAlpha,
        mixing=vdW1fRule,
        activity=nothing,
        translation=NoTranslation,
        userlocations=String[],
        ideal_userlocations=String[],
        alpha_userlocations=String[],
        activity_userlocations=String[],
        mixing_userlocations=String[],
        translation_userlocations=String[],
        verbose=false,
        assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `a`: Single Parameter (`Float64`) - Atraction Parameter
- `b`: Single Parameter (`Float64`) - Covolume
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`) - Mixed Atraction Parameter
- `b`: Pair Parameter (`Float64`) - Mixed Covolume
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Input models
- `idealmodel`: Ideal Model
- `cubicmodel`: Cubic Model

## Description

simplified CPA

## References
1. Kontogeorgis, G. M., Michelsen, M. L., Folas, G. K., Derawi, S., von Solms, N., & Stenby, E. H. (2006). Ten years with the CPA (cubic-plus-association) equation of state. Part 1. Pure compounds and self-associating systems. Industrial & Engineering Chemistry Research, 45(14), 4855–4868. [doi:10.1021/ie051305v](https://doi.org/10.1021/ie051305v)
"""
function sCPA(components;
            idealmodel=BasicIdeal,
            cubicmodel=RK,
            alpha=sCPAAlpha,
            mixing=vdW1fRule,
            activity=nothing,
            translation=NoTranslation,
            userlocations=String[],
            ideal_userlocations=String[],
            alpha_userlocations=String[],
            activity_userlocations=String[],
            mixing_userlocations=String[],
            translation_userlocations=String[],
            verbose=false,
            assoc_options = AssocOptions())

    return CPA(components;
        idealmodel=idealmodel,
        radial_dist = radial_dist,
        cubicmodel=cubicmodel,
        alpha=alpha,
        mixing=mixing,
        activity=activity,
        translation=translation,
        userlocations=userlocations,
        ideal_userlocations=ideal_userlocations,
        alpha_userlocations=alpha_userlocations,
        activity_userlocations=activity_userlocations,
        mixing_userlocations=mixing_userlocations,
        translation_userlocations=translation_userlocations,
        verbose=verbose,
        assoc_options = assoc_options)
end

export sCPA
