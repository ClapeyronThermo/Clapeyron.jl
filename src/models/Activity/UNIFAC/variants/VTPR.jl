abstract type VTPRUNIFACModel <: UNIFACModel end

struct VTPRUNIFAC{c<:EoSModel,G} <: VTPRUNIFACModel
    components::Array{String,1}
    groups::GroupParam{G}
    params::UNIFACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export VTPRUNIFAC

"""
    VTPRUNIFACModel <: UNIFACModel

    VTPRUNIFAC(components;
    puremodel = BasicIdeal,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `B`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `C`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
UNIFAC (UNIQUAC Functional-group Activity Coefficients) activity model.

Modified UNIFAC (Dortmund) implementation, only residual part, activity model used for the Volume-Translated Peng-Robinson (VTPR) EoS.

The residual part iterates over groups instead of components.

```
Gᴱ = nRT(gᴱ(res))
gᴱ(res) = -v̄∑XₖQₖlog(∑ΘₘΨₘₖ)
v̄ = ∑∑xᵢνᵢₖ for k ∈ groups,  for i ∈ components
Xₖ = (∑xᵢνᵢₖ)/v̄ for i ∈ components
Θₖ = QₖXₖ/∑QₖXₖ
Ψₖₘ = exp(-(Aₖₘ + BₖₘT + CₖₘT²)/T)
```

## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. "doi:[10.1021/i260064a004](https://doi.org/10.1021/i260064a004)"
2. Weidlich, U.; Gmehling, J. A modified UNIFAC model. 1. Prediction of VLE, hE, and.gamma..infin. Ind. Eng. Chem. Res. 1987, 26, 1372–1381.
3. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4)
"""
VTPRUNIFAC

function VTPRUNIFAC(components;
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    _components = format_gccomponents(components)
    groups = GroupParam(_components, ["Activity/UNIFAC/VTPR/VTPR_groups.csv"];group_userlocations = group_userlocations, verbose = verbose)

    params = getparams(groups, ["Activity/UNIFAC/VTPR/VTPR_like.csv", "Activity/UNIFAC/VTPR/VTPR_unlike.csv"]; userlocations = userlocations,  asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose = verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    Q  = params["Q"]
    R = deepcopy(Q)
    R.values .= 0
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String["10.1016/S0378-3812(01)00626-4"]
    model = VTPRUNIFAC(groups.components,groups,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

function excess_gibbs_free_energy(model::VTPRUNIFACModel,p,T,z)
    return excess_g_res(model,p,T,z)
end
