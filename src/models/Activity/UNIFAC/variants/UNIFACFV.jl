#a GC averaged UNIFAC.
struct UNIFACFVCache{T} <: EoSModel
    components::Vector{String}
    r::Vector{T} #molar volume
    q::Vector{T} #molar area
    qp::Vector{T} #fh volume
    Mw::Vector{T}
end

Base.eltype(::Type{UNIFACFVCache{T}}) where T = T
Base.eltype(::UNIFACFVCache{T}) where T = T

UNIFACFVCache(components,r,q,m,Mw) = UNIFACFVCache{eltype(r)}(components,r,q,m,Mw)

UNIFACFVCache(groups,params) = UNIFACFVCache(groups,params.Q,params.R,params.Mw)

function UNIFACFVCache(groups::GroupParam,Q,R,Mw)
    Mw = group_sum(groups,Mw.values)
    r = group_sum(groups,R.values)
    q = group_sum(groups,Q.values)
    qp = r .^ (3/4)
    return UNIFACFVCache(groups.components,r,q,qp,Mw)
end

function recombine_unifac_cache!(cache::UNIFACFVCache,groups,params)
    Q = params.Q
    R = params.R
    Mw = params.Mw
    group_sum!(cache.Mw,groups,Mw.values)
    group_sum!(cache.r,groups,R.values)
    #cache.r ./= Mw
    group_sum!(cache.q,groups,Q.values)
    cache.qp .= cache.r .^ (3/4)
    #cache.q ./= Mw
    return cache
end

struct UNIFACFVParam{T} <: ParametricEoSParam{T}
    volume::SingleParam{T}
    A::PairParam{T}
    R::SingleParam{T}
    Q::SingleParam{T}
    Mw::SingleParam{T}
end

UNIFACFVParam(volume,A,R,Q,Mw) = build_parametric_param(UNIFACFVParam,volume,c,A,R,Q,Mw)

abstract type UNIFACFVModel <: ActivityModel end

struct UNIFACFV{c<:EoSModel,T} <: UNIFACFVModel
    components::Array{String,1}
    groups::GroupParam{T}
    params::UNIFACFVParam{T}
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACFVCache{T}
end

function UNIFACFV(components,groups,params,puremodel,references,unifac_cache)
    c = eltype(puremodel)
    T = eltype(params)
    return UNIFACFV{c,T}(components,groups,params,puremodel,references,unifac_cache)
end

export UNIFACFV

"""
    UNIFACFVModel <: ActivityModel

    UNIFACFV(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters
- `volume`: Single Parameter (`Float64`)  - specific volume of species `[g·cm⁻³]`
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Waals volume
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asyme  trical, defaults to `0`) - Binary group Interaction Energy Parameter
- `Mw`: Single Parameter (`Float64`) - Molecular weight of groups

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties.

## Description
UNIFAC-FV (UNIFAC Free Volume) activity model. It adds a free volume term that improves the capabilities of the model to describe solvent-polymer mixtures.

The Combinatorial part corresponds to an GC-averaged modified [`UNIQUAC`](@ref) model.

```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res) + gᴱ(FV))
```
## References

1. Oishi, T., & Prausnitz, J. M. (1978). Estimation of solvent activities in polymer solutions using a group-contribution method. Industrial & Engineering Chemistry Process Design and Development, 17(3), 333–339. [doi:10.1021/i260067a021]()
"""
UNIFACFV

function UNIFACFV(components;
    puremodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    params_species = getparams(components, ["Activity/UNIFAC/UNIFACFV/UNIFACFV_like.csv"]; userlocations = userlocations, verbose = verbose, ignore_headers = ["dipprnumber","smiles","c"])

    groups = GroupParam(components, ["Activity/UNIFAC/UNIFACFV/UNIFACFV_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)
    components = groups.components
    params = getparams(groups, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_like.csv", "Activity/UNIFAC/ogUNIFAC/ogUNIFAC_unlike.csv"]; userlocations = userlocations, asymmetricparams=["A"], ignore_missing_singleparams=["A"], verbose = verbose)

    A  = params["A"]
    R  = params["R"]
    Q  = params["Q"]
    Mw = params["Mw"]
    volume  = params_species["volume"]
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = UNIFACFVParam(volume,A,R,Q,Mw)
    references = String["10.1021/i260067a021"]
    cache = UNIFACFVCache(groups,packagedparams)
    model = UNIFACFV(components,groups,packagedparams,_puremodel,references,cache)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

function recombine_impl!(model::UNIFACFVModel)
    recombine_unifac_cache!(model.unifac_cache,model.groups,model.params)
    recombine!(model.puremodel)
    return model
end

mw(model::UNIFACFVModel) = model.unifac_cache.Mw

function Ψ(model::UNIFACFVModel,V,T,z)
    A = model.params.A.values
    return @. exp(-A/T)
end

excess_g_FV(model::UNIFACFVModel,V,T,z) = excess_g_FV(model,V,T,z,FillArrays.Fill(1.1,length(model)))

function excess_g_FV(model::UNIFACFVModel,V,T,z,c)
    res = zero(Base.promote_eltype(model,V,T,z))
    Mw = model.unifac_cache.Mw
    b = 1.28
    v = model.params.volume.values
    r = model.unifac_cache.r
    ṽ = zero(res)
    r̃ = zero(res)
    for i in eachindex(z)
        zi,mi = z[i],Mw[i]
        ṽ += v[i]*zi*mi
        r̃ += r[i]*zi
    end


    #=
    The original free volume activity term proposed by Oishi and Prausnitz is not a consistent activity coefficient term (assumes that `∂v̄ₘ∂zᵢ == 0`).
    This assumption creates an inconsistent activity coefficient model (the jacobian of the activity coefficients is not symmetric.)
    To reproduce the inconsistentcy, we can mark v̄ₘ to have no derivatives, via `Clapeyron.Solvers.primalval`
    =#

    v̄ₘ = primalval(ṽ/(15.17*b*r̃))
    r̄ₘ = cbrt(v̄ₘ)
    for i in eachindex(z)
        v̄ᵢ = Mw[i]*v[i]/(15.17*b*r[i])
        r̄ᵢ,cᵢ = cbrt(v̄ᵢ),c[i]
        lnγi = 3*cᵢ*log((r̄ᵢ - 1)/(r̄ₘ - 1)) - cᵢ*((v̄ᵢ/v̄ₘ - 1)/(1 - 1/r̄ᵢ))
        res += z[i]*lnγi
    end
    return Rgas(model)*T*res
end

function excess_g_comb(model::UNIFACFVModel,V,T,z)
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    return Rgas(model)*T*gE_rt_UNIQUAC(z,r,q)
end

function excess_g_res(model::UNIFACFVModel,p,T,z)
    Ψij = Ψ(model,p,T,z)
    Q = model.params.Q.values
    return Rgas(model)*T*excess_g_res_unifac(model.groups,Q,Ψij,z)
end

function excess_gibbs_free_energy(model::UNIFACFVModel,p,T,z)
    return excess_g_comb(model,p,T,z) + excess_g_res(model,p,T,z) + excess_g_FV(model,p,T,z)
end

#old, buggy inconsistent version oof the combinatorial excess gibbs energy
function lnγ_comb_old(model::UNIFACFVModel, p, T, z)
    Mw  = model.unifac_cache.Mw
    zmw = dot(z, Mw)
    w   = z .* Mw ./ zmw
    x   = z ./ sum(z)
    r   = model.unifac_cache.r ./ Mw
    q   = model.unifac_cache.q ./ Mw
    Φ   = w .* r ./ dot(w, r)
    θ   = w .* q ./ dot(w, q)
    return @. log(Φ/x) + (1 - Φ) - 5*Mw*q*(log(Φ/θ) + (1 - Φ/θ))
end
