"""
    cPR(components::Vector{String}; idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

consistent Peng Robinson equation of state. it uses the following models:
- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

## References
1. Bell, I. H., Satyro, M., & Lemmon, E. W. (2018). Consistent twu parameters for more than 2500 pure fluids from critically evaluated experimental data. Journal of Chemical and Engineering Data, 63(7), 2402–2409. [doi:10.1021/acs.jced.7b00967](https://doi.org/10.1021/acs.jced.7b00967)
"""
function cPR(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = TwuAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation=NoTranslation,
    userlocations=String[],
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)


    #just read once if allowed.

    userlocations_cpr = String[]
    append!(userlocations_cpr,userlocations)
    if alpha === TwuAlpha
        userlocations_cpr = userlocation_merge(userlocations_cpr,userlocations)
    end

    params = getparams(components, ["cubic/cPR/cPR_single.csv"];
    userlocations=userlocations_cpr,
    verbose=verbose,
    ignore_missing_singleparams = ["v_shift","ZRA","acentricfactor"])

    k  = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]

    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",components,zeros(length(components)))
    b = PairParam("b",components,zeros(length(components)))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    if alpha !== TwuAlpha
        init_alpha = init_alphamodel(alpha,components,w,alpha_userlocations,verbose)
    else
        M = params["M"]
        N = params["N"]
        L = params["L"]
        packagedparams = TwuAlphaParam(M,N,L)
        init_alpha = TwuAlpha(packagedparams, verbose=verbose)
    end
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    packagedparams = PRParam(a,b,Tc,pc,Mw)
    references = String["10.1021/acs.jced.7b00967"]
    model = PR(components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    return model
end

export cPR