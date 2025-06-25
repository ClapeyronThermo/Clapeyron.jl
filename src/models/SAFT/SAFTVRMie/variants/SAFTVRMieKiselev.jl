"""
    SAFTVRMieKiselev::CrossOver

    SAFTVRMieKiselev(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
- `Tc`: Single Parameter (`Float64`) - critical temperature `[K]`
- `Vc`: Single Parameter (`Float64`) - critical volume `[m^3]`
- `d1`: Single Parameter (`Float64`)
- `v1`: Single Parameter (`Float64`)
- `Gi`: Single Parameter (`Float64`)
- `a20`: Single Parameter (`Float64`)
- `a21`: Single Parameter (`Float64`)


## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT-VR with Mie potential and Kiselev critical Cross-over correction.

## References
1. Lafitte, T., Apostolakou, A., Avendaño, C., Galindo, A., Adjiman, C. S., Müller, E. A., & Jackson, G. (2013). Accurate statistical associating fluid theory for chain molecules formed from Mie segments. The Journal of Chemical Physics, 139(15), 154504. [doi:10.1063/1.4819786](https://doi.org/10.1063/1.4819786)
2. Dufal, S., Lafitte, T., Haslam, A. J., Galindo, A., Clark, G. N. I., Vega, C., & Jackson, G. (2015). The A in SAFT: developing the contribution of association to the Helmholtz free energy within a Wertheim TPT1 treatment of generic Mie fluids. Molecular Physics, 113(9–10), 948–984. [doi:10.1080/00268976.2015.1029027](https://doi.org/10.1080/00268976.2015.1029027)
"""
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
    basemodel = SAFTVRMie(_components,params;idealmodel,ideal_userlocations,assoc_options,reference_state,verbose)
    return CrossOver(basemodel,critmodel;verbose)
end

default_locations(::typeof(SAFTVRMieKiselev)) = ["SAFT/SAFTVRMie/SAFTVRMieKiselev"]

export SAFTVRMieKiselev