"""
    TillnerRothFriend::MultiFluid
    TillnerRothFriend(components = ["water","ammonia"],
    Rgas =  R̄,
    verbose = false)

## Input parameters
none

##Description

Tillner-Roth and Friend model for water-ammonia mixtures.

## References
1. IAPWS G4-01 (2001). Guideline on the IAPWS Formulation 2001 for the Thermodynamic Properties of Ammonia-Water Mixtures

"""
function TillnerRothFriend(components = ["water","ammonia"],Rgas =  R̄, verbose = false)
    water = findfirst(isequal("water"),components)
    watermodel = IAPWS95()
    ammoniamodel = SingleFluid("ammonia",userlocations = ["@DB/Empiric/TLF/ammonia.json"],coolprop_userlocations = false)

    if water == 1
        pures = [watermodel,ammoniamodel]
    else
        pures = [ammoniamodel,watermodel]
    end
    specialcomp = SpecialComp(components,["ammonia"])
    mixing = TillnerRothFriendMixing(components,specialcomp)
    departure = TillnerRothFriendDeparture(components,specialcomp)
    params = MultiFluidParam(components,pures)
    references = ["IAPWS G4-01"]
    return MultiFluid(components,params,pures,mixing,departure,Rgas,references)
end

const TillnerRothModel = MultiFluid{EmpiricAncillary, TillnerRothFriendMixing, TillnerRothFriendDeparture}

export TillnerRothFriend,TillnerRothModel