"""
    QCPR(components;
        idealmodel = BasicIdeal,
        userlocations = String[], 
        ideal_userlocations = String[],
        alpha_userlocations = String[],
        mixing_userlocations = String[],
        activity_userlocations = String[],
        translation_userlocations = String[],
        reference_state = nothing,
        verbose = false)

Quantum-corrected Peng Robinson equation of state. it uses the following models:
- Translation Model: [`ConstantTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`QCPRRule`](@ref)
## References
1. Aasen, A., Hammer, M., Lasala, S., Jaubert, J.-N., & Wilhelmsen, Ø. (2020). Accurate quantum-corrected cubic equations of state for helium, neon, hydrogen, deuterium and their mixtures. Fluid Phase Equilibria, 524(112790), 112790. [doi:10.1016/j.fluid.2020.112790](https://doi.org/10.1016/j.fluid.2020.112790)
"""
function QCPR(components;
    idealmodel = BasicIdeal,
    userlocations = String[], 
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    QCPR_userlocations = userlocation_merge(["@REMOVEDEFAULTS","@DB/cubic/QCPR/QCPR_critical.csv", "@DB/cubic/QCPR/QCPR_unlike.csv"],userlocations)
    QCPR_alpha_userlocations = userlocation_merge(["@REMOVEDEFAULTS","@DB/cubic/QCPR/Twu_QCPR.csv"],alpha_userlocations)
    QCPR_translation_userlocations = userlocation_merge(["@REMOVEDEFAULTS","@DB/cubic/QCPR/QCPR_translation.csv"],translation_userlocations)

    model = PR(components;
        idealmodel = idealmodel,
        alpha = TwuAlpha,
        mixing = QCPRRule,
        activity = nothing,
        translation = ConstantTranslation,
        userlocations = QCPR_userlocations,
        ideal_userlocations = ideal_userlocations,
        alpha_userlocations = QCPR_alpha_userlocations,
        mixing_userlocations = mixing_userlocations,
        activity_userlocations = activity_userlocations,
        translation_userlocations = QCPR_translation_userlocations,
        reference_state = reference_state,
        verbose = verbose)
    setreferences!(model,String["10.1021/I160057A011"])
    return model
end

export QCPR
