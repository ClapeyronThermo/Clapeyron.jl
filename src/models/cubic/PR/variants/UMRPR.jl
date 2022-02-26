"""
    UMRPR(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = MTAlpha,
    mixing = UMRRule,
    activity = UNIFAC,
    translation=MTTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    verbose=false)

Universal Mixing Rule Peng Robinson equation of state. it uses the following models:

- Translation Model: `MTTranslation`
- Alpha Model: `MTAlpha`
- Mixing Rule Model: `UMRRule` with `UNIFAC` activity

## References

1. Voutsas, E., Magoulas, K., & Tassios, D. (2004). Universal mixing rule for cubic equations of state applicable to symmetric and asymmetric systems: Results with the Peng−Robinson equation of state. Industrial & Engineering Chemistry Research, 43(19), 6238–6246. doi:10.1021/ie049580p

"""
function UMRPR(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = MTAlpha,
    mixing = UMRRule,
    activity = UNIFAC,
    translation=MTTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    verbose=false)

    return PR(components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing=mixing,
    activity = activity,
    translation=translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    verbose = verbose)
end
export UMRPR