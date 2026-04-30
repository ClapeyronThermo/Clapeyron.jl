"""
    tcPR(components;
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

Translated and consistent Peng Robinson equation of state. It uses the following models:
- Translation Model: [`ConstantTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

If Twu parameters are not provided, they can be estimated from the acentric factor (`acentricfactor`). If translation is not provided, it can be estimated, using Rackett compresibility Factor (`ZRA`) or the acentric factor (`acentricfactor`).

## References
1. Le Guennec, Y., Privat, R., & Jaubert, J.-N. (2016). Development of the translated-consistent tc-PR and tc-RK cubic equations of state for a safe and accurate prediction of volumetric, energetic and saturation properties of pure compounds in the sub- and super-critical domains. Fluid Phase Equilibria, 429, 301–312. [doi:10.1016/j.fluid.2016.09.003](http://dx.doi.org/10.1016/j.fluid.2016.09.003)
2. Pina-Martinez, A., Le Guennec, Y., Privat, R., Jaubert, J.-N., & Mathias, P. M. (2018). Analysis of the combinations of property data that are suitable for a safe estimation of consistent twu α-function parameters: Updated parameter values for the translated-consistent tc-PR and tc-RK cubic equations of state. Journal of Chemical and Engineering Data, 63(10), 3980–3988. [doi:10.1021/acs.jced.8b00640](http://dx.doi.org/10.1021/acs.jced.8b00640)
3. Piña-Martinez, A., Privat, R., & Jaubert, J.-N. (2022). Use of 300,000 pseudo‐experimental data over 1800 pure fluids to assess the performance of four cubic equations of state: SRK , PR , tc ‐RK , and tc ‐PR. AIChE Journal. American Institute of Chemical Engineers, 68(2). [doi:10.1002/aic.17518](https://doi.org/10.1021/acs.iecr.1c03003)
"""
function tcPR(components;
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
    reference_state = nothing,
    verbose = false,
    estimate_alpha = true,
    estimate_translation = true)

    formatted_components = format_components(components)

    #just read once if allowed.
    userlocations_tcpr = String[]
    append!(userlocations_tcpr,userlocations)
    if alpha === TwuAlpha
        userlocations_tcpr = userlocation_merge(userlocations_tcpr,userlocations)
    end

    if translation == ConstantTranslation
        userlocations_tcpr = userlocation_merge(userlocations_tcpr,translation_userlocations)
    end

    params = getparams(formatted_components, ["cubic/tcPR/tcPR_single.csv"];
    userlocations = userlocations_tcpr,
    verbose = verbose,
    ignore_missing_singleparams = ["v_shift","ZRA","acentricfactor","M","N","L"])

    n = length(formatted_components)
    k  = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",formatted_components,zeros(n,n),ones(Bool,n,n))
    b = PairParam("b",formatted_components,zeros(n,n),ones(Bool,n,n))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    w = get(params,"acentricfactor",nothing)
    zra = get(params,"ZRA",nothing)

    if alpha !== TwuAlpha
        init_alpha = init_alphamodel(alpha,components,w,alpha_userlocations,verbose)
    else
        M = params["M"]
        N = params["N"]
        L = params["L"]
        for i in 1:n
            if M.ismissingvalues[i] && N.ismissingvalues[i] && L.ismissingvalues[i] && w !== nothing
                if !estimate_alpha
                    any(M.ismissingvalues) && SingleMissingError(M)
                    any(L.ismissingvalues) && SingleMissingError(L)
                    any(N.ismissingvalues) && SingleMissingError(N)
                end
                if !w.ismissingvalues[i]
                    wi = w[i]
                    #tc-PR-2020 parameters
                    Li = 0.0297*wi*wi +0.7536*wi + 0.0544
                    Mi = 0.1401*wi*wi - 0.1785*wi + 0.8678
                    L[i] = Li
                    N[i] = 2
                    M[i] = Mi
                else
                    throw(error("cannot initialize TwuAlpha in tc-PR. acentric factor not provided"))
                end
            elseif !M.ismissingvalues[i] && !N.ismissingvalues[i] && !L.ismissingvalues[i]
                #default
            else
                throw(error("invalid specification for TwuAlpha in tc-PR."))
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
        init_translation = ConstantTranslation(formatted_components,packagedparams,String[])
        #try to initialize translation if missing
        cc = init_translation.params.v_shift
        for i in 1:n
            if cc.ismissingvalues[i] 
                
                Tci = Tc[i]
                Pci = pc[i]
                R = Rgas()
                RTp = (R*Tci/Pci)
                if zra !== nothing && !zra.ismissingvalues[i]
                    cc[i] = RTp*(0.1398 - 0.5294*zra[i])
                
                elseif w !== nothing && !w.ismissingvalues[i]
                    cc[i] = RTp*(-0.0065 + 0.0198*w[i])
                else
                    throw(error("cannot initialize translation in tcPR. acentric factor or ZRA not provided"))
                end
            end
        end
    end

    packagedparams = PRParam(a,b,Tc,pc,Mw)
    references = String["10.1016/j.fluid.2016.09.003","10.1021/acs.jced.8b00640","10.1002/aic.17518","10.1021/acs.iecr.1c03003"]
    model = PR(formatted_components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
    return model
end

export tcPR