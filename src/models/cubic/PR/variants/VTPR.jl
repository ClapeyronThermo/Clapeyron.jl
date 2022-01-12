#just a function, no struct
function VTPR(components::Vector{String}; idealmodel=BasicIdeal,
    mixing = VTPRRule,
    alpha = TwuAlpha,
    translation = RackettTranslation,
    activity = VTPRUNIFAC,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
     verbose=false)

     return PR(components;
     idealmodel = idealmodel,
     alpha = alpha,
     mixing=mixing,
     activity = activity,
     translation=translation,
     ideal_userlocations = ideal_userlocations,
     alpha_userlocations = alpha_userlocations,
     mixing_userlocations = mixing_userlocations,
     verbose = verbose)
end
export VTPR