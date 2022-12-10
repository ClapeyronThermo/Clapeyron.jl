"""
    SRK(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = SoaveAlpha,
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
Soave-Redlich-Kwong equation of state. it uses the following models:

- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`SoaveAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

## References
1. Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state. Chemical Engineering Science, 27(6), 1197–1203. [doi:10.1016/0009-2509(72)80096-4](https://doi.org/10.1016/0009-2509(72)80096-4)
"""
function SRK(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = SoaveAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false, kwargs...)

     return RK(components;
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
export SRK

