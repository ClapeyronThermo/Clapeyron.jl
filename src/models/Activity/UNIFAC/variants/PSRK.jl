"""
    PSRKUNIFAC(components::Vector{String};
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Vals volume
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
    verbose = false, kwargs...)

    PRSK_userlocations = vcat("@REMOVEDEFAULTS","@DB/Activity/UNIFAC/PSRK",userlocations)
    PRSK_group_userlocations = vcat("@REMOVEDEFAULTS","@DB/Activity/UNIFAC/PSRK/PSRK_groups.csv",group_userlocations)
    model =  UNIFAC(components,
    puremodel = puremodel,
    userlocations = PRSK_userlocations,
    group_userlocations = PRSK_group_userlocations,
    pure_userlocations = pure_userlocations,
    verbose = verbose
    )
    setreferences!(model,String["10.1021/i260064a004","10.1016/j.fluid.2004.11.002"])
    return model
end

export PSRKUNIFAC
