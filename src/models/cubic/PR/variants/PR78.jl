"""
    PR78(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = PR78Alpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    verbose=false)

Peng Robinson (1978) equation of state. it uses the following models:

- Translation Model: `PR78Alpha`
- Alpha Model: `PR78Alpha`
- Mixing Rule Model: `vdW1fRule`

## References

1. Robinson DB, Peng DY. The characterization of the heptanes and heavier fractions for the GPA Peng-Robinson programs. Tulsa: Gas Processors Association; 1978
"""
function PR78(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = PR78Alpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation=NoTranslation,
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
export PR78