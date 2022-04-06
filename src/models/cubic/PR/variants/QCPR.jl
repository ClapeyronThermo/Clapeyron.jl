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

    userlocations = vcat(alpha_userlocations,translation_userlocations,userlocations)

    params = getparams(components, ["cubic/QCPR/QCPR_critical.csv", "cubic/QCPR/QCPR_unlike.csv","cubic/QCPR/Twu_QCPR.csv","cubic/QCPR/QCPR_translation.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(QCPRRule,components,nothing,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(PR,init_mixing,Tc,pc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    M = params["M"]
    N = params["N"]
    L = params["L"]
    packagedparams = TwuAlphaParam(M,N,L)
    init_alpha = TwuAlpha(packagedparams, verbose=verbose)

    c = params["c"]
    packagedparams = ConstantTranslationParam(c)
    init_translation = ConstantTranslation(packagedparams, verbose=verbose)

    icomponents = 1:length(components)
    packagedparams = PRParam(a,b,Tc,pc,Mw)
    references = String["10.1021/I160057A011"]
    model = PR(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

export QCPR
