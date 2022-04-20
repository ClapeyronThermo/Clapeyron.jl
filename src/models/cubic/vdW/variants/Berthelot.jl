"""
    Berthelot(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = ClausiusAlpha,
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

## Description
Berthelot equation of state. it uses `vdW` with the following models:

- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`ClausiusAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

## References
1. Berthelot, D. (1899). Sur une méthode purement physique pour la détermination des poids moléculaires des gaz et des poids atomiques de leurs éléments. Journal de Physique Théorique et Appliquée, 8(1), 263–274. doi:10.1051/jphystap:018990080026300
"""
function Berthelot(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = ClausiusAlpha,
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

     return vdW(components;
     idealmodel = idealmodel,
     alpha = alpha,
     mixing=mixing,
     activity = activity,
     translation = translation,
     ideal_userlocations = ideal_userlocations,
     alpha_userlocations = alpha_userlocations,
     mixing_userlocations = mixing_userlocations,
     activity_userlocations = activity_userlocations,
     translation_userlocations = translation_userlocations,
     verbose = verbose)
end
export Berthelot

