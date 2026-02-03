struct PSRKUNIFAC{c<:EoSModel,T} <: UNIFACModel
    components::Array{String,1}
    groups::GroupParam{T}
    params::UNIFACParam{T}
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache
end

default_locations(::Type{PSRKUNIFAC}) = ["Activity/UNIFAC/PSRK/PSRK_like.csv", "Activity/UNIFAC/PSRK/PSRK_unlike.csv"]


"""
    PSRKUNIFAC(components;
    puremodel = BasicIdeal,
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
UNIFAC (UNIQUAC Functional-group Activity Coefficients) activity model.

Modified [UNIFAC](@ref) (Dortmund) implementation, with parameters tuned to the Predictive Soave-Redlich-Kwong (PSRK) EoS.

## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. [doi:10.1021/i260064a004](https://doi.org/10.1021/i260064a004)
2. Weidlich, U.; Gmehling, J. A modified UNIFAC model. 1. Prediction of VLE, hE, and.gamma..infin. Ind. Eng. Chem. Res. 1987, 26, 1372–1381.
3. Horstmann, S., Jabłoniec, A., Krafczyk, J., Fischer, K., & Gmehling, J. (2005). PSRK group contribution equation of state: comprehensive revision and extension IV, including critical constants and α-function parameters for 1000 components. Fluid Phase Equilibria, 227(2), 157–164. [doi:10.1016/j.fluid.2004.11.002](https://doi.org/10.1016/j.fluid.2004.11.002)"
"""
function PSRKUNIFAC(components;
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    _components = format_gccomponents(components)
    groups = GroupParam(_components, ["Activity/UNIFAC/PSRK/PSRK_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)
    params = getparams(groups, default_locations(PSRKUNIFAC);
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
    references = String["10.1021/i260064a004","10.1016/j.fluid.2004.11.002"]
    cache = UNIFACCache(groups,packagedparams)
    model = PSRKUNIFAC(groups.components,groups,packagedparams,_puremodel,references,cache)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

#excess_g_comb(model::UNIFACModel,p,T,z=SA[1.0]) = excess_g_comb_original(model,p,T,z)

export PSRKUNIFAC
