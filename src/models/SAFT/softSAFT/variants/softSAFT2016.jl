"""
    softSAFT2016Model <: softSAFTModel

    softSAFT2016(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Description

Variant of `softSAFT` that uses the `Thol2016LJ` model instead of the default  (`Johnson1993LJ`).
"""
function softSAFT216(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

    ljmodel = Johnson1993LJ()
    return softSAFT(components;idealmodel,ljmodel,userlocations,ideal_userlocations,reference_state,verbose,assoc_options)
end

export softSAFT2016
