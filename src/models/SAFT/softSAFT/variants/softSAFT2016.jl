"""
    softSAFT2016::softSAFT

    softSAFT2016(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Description

Variant of `softSAFT` that uses the [`Thol2016LJ`](@ref) model instead of the default  ([`Thol2016LJ`](@ref)).
"""
function softSAFT2016(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

    ljmodel = Thol2016LJ()
    return softSAFT(components;idealmodel,ljmodel,userlocations,ideal_userlocations,reference_state,verbose,assoc_options)
end

export softSAFT2016
