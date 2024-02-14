function PSRKUNIFAC end
"""
    function PSRK(components;
    idealmodel = BasicIdeal,
    alpha = SoaveAlpha,
    mixing = PSRKRule,
    activity = PSRKUNIFAC,
    translation = PenelouxTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Description
Predictive Soave-Redlich-Kwong equation of state. it uses the following models:

- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`SoaveAlpha`](@ref)
- Mixing Rule Model: [`PSRKRule`](@ref) with [`PSRKUNIFAC`](@ref) activity model

##  References

1. Horstmann, S., Jabłoniec, A., Krafczyk, J., Fischer, K., & Gmehling, J. (2005). PSRK group contribution equation of state: comprehensive revision and extension IV, including critical constants and α-function parameters for 1000 components. Fluid Phase Equilibria, 227(2), 157–164. [doi:10.1016/j.fluid.2004.11.002](https://doi.org/10.1016/j.fluid.2004.11.002)
"""
function PSRK(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    activity = PSRKUNIFAC(components,
    userlocations = activity_userlocations,
    group_userlocations = group_userlocations,
    verbose = verbose)

    _components = activity.groups.components #extract pure component list

    alpha = SoaveAlpha
    mixing = PSRKRule
    translation = PenelouxTranslation

    return RK(_components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing=mixing,
    activity = activity,
    translation = translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    activity_userlocations = activity_userlocations,
    translation_userlocations = translation_userlocations,
    reference_state = reference_state,
    verbose = verbose)
end
export PSRK
