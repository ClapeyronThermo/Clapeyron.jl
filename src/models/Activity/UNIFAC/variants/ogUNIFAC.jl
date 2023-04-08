struct ogUNIFACParam <: EoSParam
    A::PairParam{Float64}
    R::SingleParam{Float64}
    Q::SingleParam{Float64}
end

abstract type ogUNIFACModel <: UNIFACModel end

struct ogUNIFAC{c<:EoSModel} <: ogUNIFACModel
    components::Array{String,1}
    groups::GroupParam
    params::ogUNIFACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache
end

@registermodel ogUNIFAC
export ogUNIFAC

"""
    ogUNIFACModel <: UNIFACModel

    ogUNIFAC(components::Vector{String};
    puremodel = PR, 
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Vals volume
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

UNIFAC (UNIQUAC Functional-group Activity Coefficients) activity model.

Original formulation.

The Combinatorial part corresponds to an GC-averaged modified UNIQUAC model. The residual part iterates over groups instead of components.

```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res))
```

Combinatorial part:
```
gᴱ(comb) = ∑[xᵢlog(Φᵢ/xᵢ) + 5qᵢxᵢlog(θᵢ/Φᵢ)]
θᵢ = qᵢxᵢ/∑qᵢxᵢ
Φᵢ = rᵢxᵢ/∑rᵢxᵢ
rᵢ = ∑Rₖνᵢₖ for k ∈ groups
qᵢ = ∑Qₖνᵢₖ for k ∈ groups
```
Residual Part:
```
gᴱ(residual) = -v̄∑XₖQₖlog(∑ΘₘΨₘₖ)
v̄ = ∑∑xᵢνᵢₖ for k ∈ groups,  for i ∈ components
Xₖ = (∑xᵢνᵢₖ)/v̄ for i ∈ components 
Θₖ = QₖXₖ/∑QₖXₖ
Ψₖₘ = exp(-(Aₖₘ/T)
```

## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. [doi:10.1021/i260064a004](https://doi.org/10.1021/i260064a004)

"""
ogUNIFAC

function ogUNIFAC(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

    groups = GroupParam(components, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_groups.csv"];group_userlocations = group_userlocations, verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_like.csv", "Activity/UNIFAC/ogUNIFAC/ogUNIFAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A"], ignore_missing_singleparams=["A"], verbose=verbose)
    A  = params["A"]
    R  = params["R"]
    Q  = params["Q"]
    
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = ogUNIFACParam(A,R,Q)
    references = String[]
    cache = UNIFACCache(groups,packagedparams)
    model = ogUNIFAC(groups.components,groups,packagedparams,_puremodel,references,cache)
    return model
end

function excess_g_comb(model::ogUNIFACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    θm = dot(q,z)*invn
    G_comp = _0
    for i ∈ @comps
        Φi = r[i]/Φm #technically xi[i]r[i]/Φm, but it gets cancelled out (log(θi/Φi))
        θi = q[i]/θm #technically xi[i]q[i]/θm, but it gets cancelled out (log(θi/Φi))
        zi = z[i]
        G_comp += zi*log(Φi) + 5*q[i]*zi*log(θi/Φi)
    end
    return G_comp
end

function Ψ(model::ogUNIFACModel,V,T,z)
    A = model.params.A.values
    return @. exp(-A/T)
end
