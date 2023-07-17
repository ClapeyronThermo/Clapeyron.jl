"""
    tcRK(components::Vector{String}; idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

translated and consistent Redlich-Kwong equation of state. it uses the following models:
- Translation Model: [`ConstantTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

If Twu parameters are not provided, they can be estimated from the acentric factor (`acentricfactor`). If translation is not provided, it can be estimated, using Rackett compresibility Factor (`ZRA`) or the acentric factor (`acentricfactor`).

## References
1. Le Guennec, Y., Privat, R., & Jaubert, J.-N. (2016). Development of the translated-consistent tc-RK and tc-RK cubic equations of state for a safe and accurate prediction of volumetric, energetic and saturation properties of pure compounds in the sub- and super-critical domains. Fluid Phase Equilibria, 429, 301–312. [doi:10.1016/j.fluid.2016.09.003](http://dx.doi.org/10.1016/j.fluid.2016.09.003)
"""
function tcRK(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = TwuAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation=ConstantTranslation,
    userlocations=String[],
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)


    #just read once if allowed.

    userlocations_tcRK = String[]
    append!(userlocations_tcRK,userlocations)
    if alpha === TwuAlpha
        userlocations_tcRK = userlocation_merge(userlocations_tcRK,userlocations)
    end

    if translation == ConstantTranslation
        userlocations_tcRK = userlocation_merge(userlocations_tcRK,translation_userlocations)
    end

    params = getparams(components, ["cubic/tcRK/tcRK_single.csv"];
    userlocations=userlocations_tcRK,
    verbose=verbose,
    ignore_missing_singleparams = ["v_shift","ZRA","acentricfactor","M","N","L"])

    n = length(components)
    k  = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",components,zeros(length(components)))
    b = PairParam("b",components,zeros(length(components)))
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
                if !w.ismissingvalues[i]
                    wi = w[i]
                    Li = 0.0947*wi*wi + 0.6871*wi + 0.1508
                    Mi =  0.1615*wi*wi - 0.2349*wi + 0.8876
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
        init_alpha = TwuAlpha(packagedparams, verbose=verbose)
    end

    if translation !== ConstantTranslation
        init_translation = init_model(translation,components,translation_userlocations,verbose)
    else
        c = params["v_shift"]
        packagedparams = ConstantTranslationParam(c)
        init_translation = ConstantTranslation(packagedparams, verbose=verbose)
        #try to initialize translation if missing
        cc = init_translation.params.v_shift
        for i in 1:n
            if cc.ismissingvalues[i]
                Tci = Tc[i]
                Pci = Pc[i]
                R = Rgas()
                RTp = (R*Tci/Pci)
                if !zra.ismissingvalues[i]
                    cc[i] = RTp*(0.1487 - 0.5052*zra[i])
                elseif !w.ismissingvalues[i]
                    cc[i] = RTp*(0.0096 + 0.0172*w[i])
                else
                    throw(error("cannot initialize translation in tcRK. acentric factor or ZRA not provided"))
                end
            end
        end
    end

    packagedparams = RKParam(a,b,Tc,pc,Mw)
    references = String["10.1016/j.fluid.2016.09.003"]
    model = RK(components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    return model
end