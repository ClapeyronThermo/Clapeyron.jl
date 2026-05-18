struct UNIFAC2{c<:EoSModel,T} <: UNIFACModel
    components::Array{String,1}
    groups::GroupParam{T}
    params::UNIFACParam{T}
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache{T}
end

function UNIFAC2(components,groups,params,puremodel,references,unifac_cache)
    c = eltype(puremodel)
    T = eltype(params)
    return UNIFAC2{c,T}(components,groups,params,puremodel,references,unifac_cache)
end

default_locations(::Type{UNIFAC2}) = ["Activity/UNIFAC/UNIFAC_like.csv", "Activity/UNIFAC/UNIFAC2_unlike.csv"]
const modUNIFAC2 = UNIFAC2
const UNIFAC2_0 = UNIFAC2
export UNIFAC2

"""
    UNIFACModel <: ActivityModel

    UNIFAC2(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Waals volume
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `B`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `C`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
UNIFAC 2.0 activity model.
Modified UNIFAC 2.0 (Dortmund) implementation.
The method is identical to [`UNIFAC`](@ref) but with a new parameters fitted by matrix completion methods.

## References
1. Hayer, N., Hasse, H., Jirasek, F.: Modified UNIFAC 2.0-A Group-Contribution Method Completed with Machine Learning, Ind. Eng. Chem. Res. 64 (2025) 10304–10313, DOI: [10.1021/acs.iecr.5c00077](https://doi.org/10.1021/acs.iecr.5c00077).
"""
UNIFAC2

function UNIFAC2(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    groups = GroupParam(components, ["Activity/UNIFAC/UNIFAC_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)

    params = getparams(groups, default_locations(UNIFAC2);
                        userlocations = userlocations,
                        asymmetricparams=["A","B","C"],
                        ignore_missing_singleparams=["A","B","C"],
                        verbose = verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    R  = params["R"]
    Q  = params["Q"]
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String["10.48550/arXiv.2412.12962"]
    cache = UNIFACCache(groups,packagedparams)
    model = UNIFAC2(groups.components,groups,packagedparams,_puremodel,references,cache)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end
