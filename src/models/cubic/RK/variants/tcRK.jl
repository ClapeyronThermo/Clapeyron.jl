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
    verbose = false,
    estimate_alpha = true,
    estimate_translation = true)

translated and consistent Redlich-Kwong equation of state. it uses the following models:
- Translation Model: [`ConstantTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

If Twu parameters are not provided, they can be estimated from the acentric factor (`acentricfactor`). If translation is not provided, it can be estimated, using Rackett compresibility Factor (`ZRA`) or the acentric factor (`acentricfactor`).

The use of estimates for the alpha function and volume translation can be turned off by passing `estimate_alpha = false` or `estimate_translation = false`.

## References
1. Le Guennec, Y., Privat, R., & Jaubert, J.-N. (2016). Development of the translated-consistent tc-PR and tc-RK cubic equations of state for a safe and accurate prediction of volumetric, energetic and saturation properties of pure compounds in the sub- and super-critical domains. Fluid Phase Equilibria, 429, 301–312. [doi:10.1016/j.fluid.2016.09.003](http://dx.doi.org/10.1016/j.fluid.2016.09.003)
2. Pina-Martinez, A., Le Guennec, Y., Privat, R., Jaubert, J.-N., & Mathias, P. M. (2018). Analysis of the combinations of property data that are suitable for a safe estimation of consistent twu α-function parameters: Updated parameter values for the translated-consistent tc-PR and tc-RK cubic equations of state. Journal of Chemical and Engineering Data, 63(10), 3980–3988. [doi:10.1021/acs.jced.8b00640](http://dx.doi.org/10.1021/acs.jced.8b00640)
3. Piña-Martinez, A., Privat, R., & Jaubert, J.-N. (2022). Use of 300,000 pseudo‐experimental data over 1800 pure fluids to assess the performance of four cubic equations of state: SRK , PR , tc ‐RK , and tc ‐PR. AIChE Journal. American Institute of Chemical Engineers, 68(2). [doi:10.1002/aic.17518](https://doi.org/10.1021/acs.iecr.1c03003)
"""
function tcRK(components;
    idealmodel = BasicIdeal,
    alpha = TwuAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = ConstantTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    estimate_alpha = true,
    estimate_translation = true,
    reference_state = nothing,
    verbose = false)


    #just read once if allowed.

    userlocations_tcRK = String[]
    append!(userlocations_tcRK,userlocations)
    if alpha === TwuAlpha
        userlocations_tcRK = userlocation_merge(userlocations_tcRK,userlocations)
    end

    if translation == ConstantTranslation
        userlocations_tcRK = userlocation_merge(userlocations_tcRK,translation_userlocations)
    end
    formatted_components = format_components(components)
    params = getparams(formatted_components, ["cubic/tcRK/tcRK_single.csv"];
    userlocations = userlocations_tcRK,
    verbose = verbose,
    ignore_missing_singleparams = ["v_shift","ZRA","acentricfactor","M","N","L"])

    n = length(formatted_components)
    k  = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",components,zeros(n))
    b = PairParam("b",components,zeros(n))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    w = params["acentricfactor"]
    zra = params["ZRA"]

    if alpha !== TwuAlpha
        init_alpha = init_alphamodel(alpha,components,w,alpha_userlocations,verbose)
    else
        M = params["M"]
        N = params["N"]
        L = params["L"]
        for i in 1:n
            if M.ismissingvalues[i] && N.ismissingvalues[i] && L.ismissingvalues[i]
                if !estimate_alpha
                    any(M.ismissingvalues) && SingleMissingError(M)
                    any(L.ismissingvalues) && SingleMissingError(L)
                    any(N.ismissingvalues) && SingleMissingError(N)
                end
                if !w.ismissingvalues[i]
                    wi = w[i]
                    Li = 0.0611*wi*wi + 0.7535*wi + 0.1359
                    Mi = 0.1709*wi*wi - 0.2063*wi + 0.8787
                    L[i] = Li
                    N[i] = 2
                    M[i] = Mi
                else
                    throw(error("cannot initialize TwuAlpha in tc-RK. acentric factor not provided"))
                end
            elseif !M.ismissingvalues[i] && !N.ismissingvalues[i] && !L.ismissingvalues[i]
                #default
            else
                throw(error("invalid specification for TwuAlpha in tc-RK."))
            end
        end
        packagedparams = TwuAlphaParam(M,N,L)
        init_alpha = TwuAlpha(formatted_components,packagedparams,default_references(TwuAlpha))
    end

    if translation !== ConstantTranslation
        init_translation = init_model(translation,components,translation_userlocations,verbose)
    else
        c = params["v_shift"]
        packagedparams = ConstantTranslationParam(c)
        init_translation = ConstantTranslation(formatted_components, packagedparams, String[])
        #try to initialize translation if missing
        cc = init_translation.params.v_shift
        for i in 1:n
            if cc.ismissingvalues[i]
                Tci = Tc[i]
                Pci = pc[i]
                R = Rgas()
                RTp = (R*Tci/Pci)
                !estimate_translation && SingleMissingError(c)
                if !zra.ismissingvalues[i]
                    cc[i] = RTp*(0.2150 - 0.7314*zra[i])
                elseif !w.ismissingvalues[i]
                    cc[i] = RTp*(0.0227 + 0.0093*w[i])
                else
                    throw(error("cannot initialize translation in tcRK. acentric factor or ZRA not provided"))
                end
            end
        end
    end

    packagedparams = RKParam(a,b,Tc,pc,Mw)
    references = String["10.1016/j.fluid.2016.09.003","10.1021/acs.jced.8b00640","10.1002/aic.17518","10.1021/acs.iecr.1c03003"]
    model = RK(formatted_components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
    return model
end

export tcRK