function gcPCPSAFT(components;
    mixing = :hetero,
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    verbose = false,
    reference_state = nothing,
    assoc_options = AssocOptions())

    if mixing == :homo
        return HomogcPCPSAFT(components;
        idealmodel = idealmodel,
        userlocations = userlocations,
        group_userlocations = group_userlocations,
        ideal_userlocations = ideal_userlocations,
        verbose = verbose,
        reference_state = reference_state,
        assoc_options = assoc_options)
    elseif mixing == :hetero
        return HeterogcPCPSAFT(components;
        idealmodel = idealmodel,
        userlocations = userlocations,
        group_userlocations = group_userlocations,
        ideal_userlocations = ideal_userlocations,
        verbose = verbose,
        reference_state = reference_state,
        assoc_options = assoc_options)
    else
        throw(ArgumentError("gcPCSAFT: mixing must be :homo or :hetero"))
    end
end

function __gcpcpsaft_combine(mixing)
    if mixing == :homo
        return :cr1
    elseif mixing ==:hetero
        return :cr1
    else
        throw(ArgumentError("gcPCSAFT: mixing must be :homo or :hetero"))
    end
end

function gcPCSAFT(components;
    mixing = :hetero,
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    verbose = false,
    reference_state = nothing,
    assoc_options = AssocOptions(combining = __gcpcpsaft_combine(mixing)))

    if mixing == :homo
        return HomogcPCPSAFT(components;
        idealmodel = idealmodel,
        userlocations = userlocations,
        group_userlocations = group_userlocations,
        ideal_userlocations = ideal_userlocations,
        verbose = verbose,
        reference_state = reference_state,
        assoc_options = assoc_options)
    elseif mixing == :hetero
        return HeterogcPCPSAFT(components;
        idealmodel = idealmodel,
        userlocations = userlocations,
        group_userlocations = group_userlocations,
        ideal_userlocations = ideal_userlocations,
        verbose = verbose,
        reference_state = reference_state,
        assoc_options = assoc_options)
    else
        throw(ArgumentError("gcPCSAFT: mixing must be :homo or :hetero"))
    end
end

Base.@deprecate_binding gcPPCSAFT gcPCPSAFT
