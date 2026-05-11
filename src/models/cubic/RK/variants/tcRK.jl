"""
    tcRK(components::Vector{String};
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

Translated and consistent Redlich-Kwong equation of state. It uses the following models:
- Translation Model: [`ConstantTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

If Twu parameters are not provided, they can be estimated from the acentric factor (`acentricfactor`). If translation is not provided, it can be estimated, using Rackett compresibility Factor (`ZRA`) or the acentric factor (`acentricfactor`).

## References
1. Le Guennec, Y., Privat, R., & Jaubert, J.-N. (2016). Development of the translated-consistent tc-PR and tc-RK cubic equations of state for a safe and accurate prediction of volumetric, energetic and saturation properties of pure compounds in the sub- and super-critical domains. Fluid Phase Equilibria, 429, 301–312. [doi:10.1016/j.fluid.2016.09.003](http://dx.doi.org/10.1016/j.fluid.2016.09.003)
2. Pina-Martinez, A., Le Guennec, Y., Privat, R., Jaubert, J.-N., & Mathias, P. M. (2018). Analysis of the combinations of property data that are suitable for a safe estimation of consistent twu α-function parameters: Updated parameter values for the translated-consistent tc-PR and tc-RK cubic equations of state. Journal of Chemical and Engineering Data, 63(10), 3980–3988. [doi:10.1021/acs.jced.8b00640](http://dx.doi.org/10.1021/acs.jced.8b00640)
3. Piña-Martinez, A., Privat, R., & Jaubert, J.-N. (2022). Use of 300,000 pseudo‐experimental data over 1800 pure fluids to assess the performance of four cubic equations of state: SRK , PR , tc ‐RK , and tc ‐PR. AIChE Journal. American Institute of Chemical Engineers, 68(2). [doi:10.1002/aic.17518](https://doi.org/10.1021/acs.iecr.1c03003)
"""
function tcRK(components;
    idealmodel = BasicIdeal,
    alpha = tcTwuAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = tcTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    #just read once if allowed.

    formatted_components = format_components(components)
    tc_userlocations = String[]
    userlocation_merge(tc_userlocations,userlocations)

    alpha === tcTwuAlpha && userlocation_merge(tc_userlocations,userlocations)
    translation == tcTranslation && userlocation_merge(tc_userlocations,translation_userlocations)

    params = getparams(formatted_components, ["cubic/tcRK/tcRK_single.csv"];
                        userlocations = tc_userlocations,
                        verbose = verbose,
                        ignore_missing_singleparams = ["v_shift","ZRA","acentricfactor","M","N","L"])


    _alpha = if alpha == tcTwuAlpha
        alphaparams = transform_params(tcTwuAlpha,params,formatted_components)
        alpha_packagedparams = build_eosparam(tcTwuAlphaParam,alphaparams)
        init_alpha = tcTwuAlpha(formatted_components,alpha_packagedparams,default_references(tcTwuAlpha))
    else
        alpha
    end

    _translation = if translation == tcTranslation
        translationparams = transform_params(tcTranslation,params,formatted_components)
        translation_packagedparams = build_eosparam(tcTranslationParam,translationparams)
        init_alpha = tcTranslation(formatted_components,translation_packagedparams,default_references(tcTranslation))
    else
        translation
    end

    model = CubicModel(RK,params,formatted_components;
                    alpha = _alpha,translation = _translation,idealmodel,mixing,activity,
                    userlocations,ideal_userlocations,alpha_userlocations,activity_userlocations,mixing_userlocations,translation_userlocations,
                    reference_state, verbose)

    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
    references = ["10.1016/j.fluid.2016.09.003","10.1021/acs.jced.8b00640","10.1002/aic.17518","10.1021/acs.iecr.1c03003"]
    setreferences!(model,references)
    return model
end

export tcRK