abstract type VTPRUNIFACModel <: UNIFACModel end


struct VTPRUNIFACCache <: EoSModel
    components::Vector{String}
    m::Vector{Float64}
end

struct VTPRUNIFAC{c<:EoSModel} <: VTPRUNIFACModel
    components::Array{String,1}
    groups::GroupParam
    params::UNIFACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::VTPRUNIFACCache
end

@registermodel VTPRUNIFAC
export VTPRUNIFAC

"""
    VTPRUNIFACModel <: UNIFACModel

    VTPRUNIFAC(components::Vector{String};
    puremodel = BasicIdeal,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

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
    verbose = false, kwargs...)

    groups = GroupParam(components, ["Activity/UNIFAC/VTPR/VTPR_groups.csv"];group_userlocations = group_userlocations, verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/VTPR/VTPR_like.csv", "Activity/UNIFAC/VTPR/VTPR_unlike.csv"]; userlocations=userlocations,  asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    Q  = params["Q"]
    R = deepcopy(Q)
    R.values .= 0
    cache = VTPRUNIFACCache(groups)
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String["10.1016/S0378-3812(01)00626-4"]
    model = VTPRUNIFAC(groups.components,groups,packagedparams,_puremodel,references,cache)
    return model
end

function recombine_unifac_cache!(cache::VTPRUNIFACCache,groups,params)
    group_sum!(cache.m,groups,nothing)
    return cache
end

function VTPRUNIFACCache(groups::GroupParam)
    m = group_sum(groups,nothing)
    return VTPRUNIFACCache(groups.components,m)
end

function excess_gibbs_free_energy(model::VTPRUNIFACModel,V,T,z)
    lnγ = lnγ_res(model,V,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end
