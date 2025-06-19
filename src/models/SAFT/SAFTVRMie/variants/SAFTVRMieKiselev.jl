function SAFTVRMieKiselev(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    assoc_options = Clapeyron.default_assoc_options(SAFTVRMie),
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params = getparams(_components, default_locations(SAFTVRMieKiselev); userlocations = userlocations, verbose = verbose)

    critparams = build_eosparam(Kiselev2000Param,params)
    critmodel = Kiselev2000(_components,critparams,default_references(Kiselev2000))

    #build SAFTVRMie
    
    basemodel = SAFTVRMie15(_components,params;idealmodel,ideal_userlocations,assoc_options,reference_state,verbose)
    return CrossOver(basemodel,critmodel;verbose)
end

default_locations(::typeof(SAFTVRMieKiselev)) = ["SAFT/SAFTVRMie/SAFTVRMieKiselev"]

export SAFTVRMieKiselev