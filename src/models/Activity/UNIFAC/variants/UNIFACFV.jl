#a GC averaged UNIFAC.
struct UNIFACFVCache{T} <: EoSModel
    components::Vector{String}
    r::Vector{T}
    q::Vector{T}
    m::Vector{T}
    Mw::Vector{T}
end

Base.eltype(::Type{UNIFACFVCache{T}}) where T = T
Base.eltype(::UNIFACFVCache{T}) where T = T


UNIFACFVCache(components,r,q,m,Mw) = UNIFACFVCache{eltype(r)}(components,r,q,m,Mw)

UNIFACFVCache(groups,params) = UNIFACFVCache(groups,params.Q,params.R,params.Mw)

function UNIFACFVCache(groups::GroupParam,Q,R,Mw)
    Mw = group_sum(groups,Mw.values)
    r = group_sum(groups,R.values) ./ Mw
    q = group_sum(groups,Q.values) ./ Mw
    m = group_sum(groups,nothing)
    return UNIFACFVCache(groups.components,r,q,m,Mw)
end

function recombine_unifac_cache!(cache::UNIFACFVCache,groups,params)
    Q = params.Q
    R = params.R
    Mw = params.Mw
    group_sum!(cache.r,groups,R.values)
    cache.r ./= Mw
    group_sum!(cache.q,groups,Q.values)
    cache.q ./= Mw
    group_sum!(cache.m,groups,nothing)
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
    UNIFACFV_cache::UNIFACFVCache{T}
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
- `volume`: Single Parameter (`Float64`)  - specific volume of species `[g/cm^3]`
- `R`: Single Parameter (`Float64`)  - Normalized group Van der Waals volume
- `Q`: Single Parameter (`Float64`) - Normalized group Surface Area
- `A`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Binary group Interaction Energy Parameter
- `Mw`: Single Parameter (`Float64`) - Molecular weight of groups

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
UNIFACFV (UNIFAC Free Volume) activity model. specialized for solvent-polymer mixtures

The Combinatorial part corresponds to an GC-averaged modified [`UNIQUAC`](@ref) model.

```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res) + gᴱ(FV))
```
## References

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
    references = String["10.1021/i260064a004"]
    cache = UNIFACFVCache(groups,packagedparams)
    model = UNIFACFV(components,groups,packagedparams,_puremodel,references,cache)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

function recombine_impl!(model::UNIFACFVModel)
    recombine_unifac_cache!(model.UNIFACFV_cache,model.groups,model.params)
    recombine!(model.puremodel)
    return model
end

function lnγ_impl!(res,model::UNIFACFVModel,V,T,z)
    _data = @f(data)
    res .= 0
    res .+= @f(lnγ_comb,_data)
    res .+= @f(lnγ_res,_data)
    res .+= @f(lnγ_FV,_data)
    return res
end

function data(model::UNIFACFVModel,V,T,z)
    Mw = model.UNIFACFV_cache.Mw
    zmw = dot(z,Mw)
    x = z ./ sum(z)
    w = z .* Mw / zmw
    c = FillArrays.Fill(1.1,length(model))
    return w,x,c
end

function lnγ_comb(model::UNIFACFVModel,V,T,z,_data=@f(data))
    w,x = _data
    Mw = model.UNIFACFV_cache.Mw
    r =model.UNIFACFV_cache.r
    q =model.UNIFACFV_cache.q
    Φ = w.*r/dot(w,r)
    θ = w.*q/dot(w,q)
    lnγ_comb = @. log(Φ/x)+(1-Φ)-5*Mw*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end

function lnγ_res(model::UNIFACFVModel,V,T,z,_data=@f(data))
    v  = model.groups.n_flattenedgroups
    _Ψ = @f(Ψ)
    lnΓ_ = @f(lnΓ,_Ψ)
    lnΓi_ = @f(lnΓi,_Ψ)
    lnγ_res_ = [sum(v[i][k].*(lnΓ_[k].-lnΓi_[i][k]) for k ∈ @groups) for i ∈ @comps]
    return lnγ_res_
end

function lnΓ(model::UNIFACFVModel,V,T,z,_Ψ = @f(Ψ))
    Mw = model.params.Mw.values
    Q = model.params.Q.values ./ Mw
    v  = model.groups.n_flattenedgroups
    x = z ./ sum(z)
    W = sum(v[i][:].*Mw*x[i] for i ∈ @comps) ./ sum(sum(v[i][k]*Mw[k]*x[i] for k ∈ @groups) for i ∈ @comps)
    θ = W.*Q / dot(W,Q)
    lnΓ_ = Mw.*Q.*(1 .-log.(sum(θ[m]*_Ψ[m,:] for m ∈ @groups)) .- sum(θ[m]*_Ψ[:,m]./sum(θ[n]*_Ψ[n,m] for n ∈ @groups) for m ∈ @groups))
    return lnΓ_
end

function lnΓi(model::UNIFACFVModel,V,T,z,_Ψ = @f(Ψ))
    Mw = model.params.Mw.values
    Q = model.params.Q.values ./ Mw
    v  = model.groups.n_flattenedgroups
    W = [v[i][:].*Mw ./ sum(v[i][k]*Mw[k] for k ∈ @groups) for i ∈ @comps]
    θ = [W[i][:].*Q ./ sum(W[i][n]*Q[n] for n ∈ @groups) for i ∈ @comps]
    lnΓi_ = [Mw.*Q.*(1 .-log.(sum(θ[i][m]*_Ψ[m,:] for m ∈ @groups)) .- sum(θ[i][m]*_Ψ[:,m]./sum(θ[i][n]*_Ψ[n,m] for n ∈ @groups) for m ∈ @groups)) for i ∈ @comps]
    return lnΓi_
end

function Ψ(model::UNIFACFVModel,V,T,z)
    A = model.params.A.values
    return @. exp(-A/T)
end

function lnγ_FV(model::UNIFACFVModel,V,T,z,_data=@f(data))
    w,x,c = _data
    b = 1.28
    v = model.params.volume.values
    r = model.UNIFACFV_cache.r
    v̄  = @. v/(15.17*b*r)
    v̄ₘ = dot(v,w)/(15.17*b*dot(r,w))
    return @. 3*c*log((v̄^(1/3)-1)/(cbrt(v̄ₘ)-1))-c*((v̄/v̄ₘ-1)/(1-v̄^(-1/3)))
end

function excess_g_SG(model::UNIFACFVModel,p,T,z)
    lnγ = lnγ_SG(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function excess_g_res(model::UNIFACFVModel,p,T,z)
    lnγ = lnγ_res(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end