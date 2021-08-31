#just a function, no struct
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
    verbose=false)

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