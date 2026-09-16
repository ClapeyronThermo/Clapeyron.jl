struct UNIFACFVPolyParam{T} <: ParametricEoSParam{T}
    volume::SingleParam{T}
    c::SingleParam{T}
    A::PairParam{T}
    R::SingleParam{T}
    Q::SingleParam{T}
    Mw::SingleParam{T}
end

UNIFACFVPolyParam(volume,c,A,R,Q,Mw) = build_parametric_param(UNIFACFVPolyParam,volume,c,A,R,Q,Mw)

abstract type UNIFACFVPolyModel <: UNIFACFVModel end

struct UNIFACFVPoly{c<:EoSModel,T} <: UNIFACFVPolyModel
    components::Array{String,1}
    groups::GroupParam{T}
    params::UNIFACFVPolyParam{T}
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACFVCache{T}
end

function UNIFACFVPoly(components,groups,params,puremodel,references,unifac_cache)
    c = eltype(puremodel)
    T = eltype(params)
    return UNIFACFVPoly{c,T}(components,groups,params,puremodel,references,unifac_cache)
end

export UNIFACFVPoly

"""
    UNIFACFVPolyModel <: ActivityModel

    UNIFACFVPoly(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters
- `volume`: Single Parameter (`Float64`)  - specific volume of species `[g·cm⁻³]`
- `c` Single Parameter  (`Float64`)  - number of external degrees of freedom per solvent molecule
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Waals volume
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `Mw`: Single Parameter (`Float64`) - Molecular weight of groups

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
UNIFAC-FV (polymer) (UNIFAC Free Volume) activity model. Specialized for polymer blends.

The Combinatorial part corresponds to an GC-averaged Flory-Huggings contribution.

```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res) + gᴱ(FV))
```
"""
UNIFACFVPoly

function UNIFACFVPoly(components;
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    params_species = getparams(components, ["Activity/UNIFAC/UNIFACFV/UNIFACFV_like.csv"]; userlocations = userlocations, verbose = verbose)

    groups = GroupParam(components, ["Activity/UNIFAC/UNIFACFV/UNIFACFV_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)
    components = groups.components
    params = getparams(groups, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_like.csv", "Activity/UNIFAC/ogUNIFAC/ogUNIFAC_unlike.csv"]; userlocations = userlocations, asymmetricparams=["A"], ignore_missing_singleparams=["A"], verbose = verbose)
    A  = params["A"]
    R  = params["R"]
    Q  = params["Q"]
    Mw = params["Mw"]
    volume  = params_species["volume"]
    c  = params_species["c"]
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = UNIFACFVPolyParam(volume,c,A,R,Q,Mw)
    references = String["10.1021/i260064a004"]
    cache = UNIFACFVCache(groups,packagedparams)
    model = UNIFACFVPoly(components,groups,packagedparams,_puremodel,references,cache)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

excess_g_FV(model::UNIFACFVPolyModel,V,T,z) = excess_g_FV(model,V,T,z,model.params.c.values)

function excess_g_comb(model::UNIFACFVPolyModel,V,T,z)
    qp = model.unifac_cache.qp
    return Rgas(model)*T*gE_rt_FH(z,qp)
end

#old, buggy inconsistent version of the combinatorial excess gibbs energy for UNIFACFVPoly
function lnγ_comb_old(model::UNIFACFVPolyModel,V,T,z)
    Mw  = model.unifac_cache.Mw
    zmw = dot(z, Mw)
    w   = z .* Mw ./ zmw
    x   = z ./ sum(z)
    r   = model.unifac_cache.r ./ Mw
    wr = @sum(w[i]*r[i]^(3/4))
    Φ = w.*r.^(3/4)/wr
    lnγ_comb = @. log(Φ/x)+(1-Φ)
    return lnγ_comb
end