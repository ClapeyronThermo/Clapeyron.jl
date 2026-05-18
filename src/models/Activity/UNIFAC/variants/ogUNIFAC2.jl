struct ogUNIFAC2{c<:EoSModel,T} <: ogUNIFACModel
    components::Array{String,1}
    groups::GroupParam{T}
    params::ogUNIFACParam{T}
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache{T}
end

function ogUNIFAC2(components,groups,params,puremodel,references,unifac_cache)
    c = eltype(puremodel)
    T = eltype(params)
    return ogUNIFAC2{c,T}(components,groups,params,puremodel,references,unifac_cache)
end

const ogUNIFAC2_0 = ogUNIFAC2
export ogUNIFAC2

"""
    ogUNIFACModel <: UNIFACModel

    ogUNIFAC2(components;
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

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

UNIFAC 2.0 (UNIQUAC Functional-group Activity Coefficients) activity model.
Original formulation.
The method is identical to [`ogUNIFAC`](@ref) but with a new parameters fitted by matrix completion methods.

## References
1. Hayer, N., Wendel, T., Mandt, S., Hasse, H., Jirasek, F., Advancing Thermodynamic Group-Contribution Methods by Machine Learning: UNIFAC 2.0, Chemical Engineering Journal 504 (2025) 158667. [10.1016/j.cej.2024.158667](https://doi.org/10.1016/j.cej.2024.158667).

"""
ogUNIFAC2

default_locations(::Type{ogUNIFAC2}) = ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_like.csv", "Activity/UNIFAC/ogUNIFAC/ogUNIFAC2_unlike.csv"]

function ogUNIFAC2(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    groups = GroupParam(components, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)

    params = getparams(groups, default_locations(ogUNIFAC2); userlocations = userlocations, asymmetricparams=["A"], ignore_missing_singleparams=["A"], verbose = verbose)
    A  = params["A"]
    R  = params["R"]
    Q  = params["Q"]
    
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = ogUNIFACParam(A,R,Q)
    references = String["10.1016/j.cej.2024.158667"]
    cache = UNIFACCache(groups,packagedparams)
    model = ogUNIFAC2(groups.components,groups,packagedparams,_puremodel,references,cache)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end
