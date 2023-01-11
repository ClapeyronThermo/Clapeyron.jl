"""
    VTPR(components; idealmodel=BasicIdeal,
    userlocations=String[],
    group_userlocations = String[]
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)
Volume-translated Peng Robinson equation of state. it uses the following models:
- Translation Model: [`RackettTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`VTPRRule`](@ref) with [`VTPRUNIFAC`](@ref) activity
## References
1. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4)
"""
function VTPR(components;
    idealmodel=BasicIdeal,
    userlocations=String[], 
    group_userlocations = String[],
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

    activity = VTPRUNIFAC(components,
            userlocations = activity_userlocations,
            group_userlocations = group_userlocations,
            verbose = verbose)

    _components = activity.groups.components #extract pure component list

    translation = RackettTranslation
    alpha = TwuAlpha
    mixing = VTPRRule

    return PR(_components;
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
export VTPR
