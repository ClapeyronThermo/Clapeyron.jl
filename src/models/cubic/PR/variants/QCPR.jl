"""
    QCPR(components::Vector{String}; idealmodel=BasicIdeal,
        userlocations=String[], 
        ideal_userlocations=String[],
        alpha_userlocations = String[],
        mixing_userlocations = String[],
        activity_userlocations = String[],
        translation_userlocations = String[],
        verbose=false)

Quantum-corrected Peng Robinson equation of state. it uses the following models:

- Translation Model: `ConstantTranslation`
- Alpha Model: `TwuAlpha`
- Mixing Rule Model: `QCPRRule`

## References

1. Aasen, A., Hammer, M., Lasala, S., Jaubert, J.-N., & Wilhelmsen, Ø. (2020). Accurate quantum-corrected cubic equations of state for helium, neon, hydrogen, deuterium and their mixtures. Fluid Phase Equilibria, 524(112790), 112790. doi:10.1016/j.fluid.2020.112790

"""
function QCPR(components::Vector{String}; idealmodel=BasicIdeal,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

    crit = "cubic/QCPR/QCPR_critical.csv"

    twu_alpha = "cubic/QCPR/Twu_QCPR.csv"
    const_translation = "cubic/QCPR/QCPR_translation.csv"

    userlocations = vcat(crit,userlocations)
    alpha_userlocations = vcat(twu_alpha,alpha_userlocations)
   
    translation_userlocations = vcat(const_translation,translation_userlocations)

    model = PR(components; idealmodel=idealmodel,
    alpha = TwuAlpha,
    mixing = QCPRRule,
    activity=nothing,
    translation=ConstantTranslation,
    userlocations=userlocations, 
    ideal_userlocations=ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    activity_userlocations = activity_userlocations,
    translation_userlocations = translation_userlocations,
    verbose=verbose)

    push!(model.references,"10.1016/j.fluid.2020.112790")
    return model
end

export QCPR
