function CPC_RKE(components::Vector{String}; 
    idealmodel=BasicIdeal,
    cubicmodel = RK,
    alpha = RKAlpha,
    mixing = CPCRule,
    rdf = ElliottRDF,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    rdf_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

    return CubicPlusChain(components; 
    idealmodel,
    cubicmodel,
    alpha,
    mixing,
    rdf,
    translation,
    userlocations, 
    ideal_userlocations,
    alpha_userlocations,
    mixing_userlocations,
    rdf_userlocations,
    translation_userlocations,
    verbose)
end 
