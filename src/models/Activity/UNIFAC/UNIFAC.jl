struct UNIFACParam <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
    C::PairParam{Float64}
    R::SingleParam{Float64}
    Q::SingleParam{Float64}
end

abstract type UNIFACModel <: ActivityModel end

struct UNIFAC{c<:EoSModel} <: UNIFACModel
    components::Array{String,1}
    groups::GroupParam
    params::UNIFACParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache
end

@registermodel UNIFAC
const modUNIFAC = UNIFAC
export UNIFAC

"""
    UNIFACModel <: ActivityModel

    UNIFAC(components;
    puremodel = PR,
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
Modified UNIFAC (Dortmund) implementation.
The Combinatorial part corresponds to an GC-averaged modified [`UNIQUAC`](@ref) model. The residual part iterates over groups instead of components.
```
Gᴱ = nRT(gᴱ(comb) + gᴱ(res))
```
Combinatorial part:
```
gᴱ(comb) = ∑[xᵢlog(Φ'ᵢ) + 5qᵢxᵢlog(θᵢ/Φᵢ)]
θᵢ = qᵢxᵢ/∑qᵢxᵢ
Φᵢ = rᵢxᵢ/∑rᵢxᵢ
Φ'ᵢ = rᵢ^(0.75)/∑xᵢrᵢ^(0.75)
rᵢ = ∑Rₖνᵢₖ for k ∈ groups
qᵢ = ∑Qₖνᵢₖ for k ∈ groups
```
Residual Part:
```
gᴱ(residual) = -v̄∑XₖQₖlog(∑ΘₘΨₘₖ)
v̄ = ∑∑xᵢνᵢₖ for k ∈ groups,  for i ∈ components
Xₖ = (∑xᵢνᵢₖ)/v̄ for i ∈ components
Θₖ = QₖXₖ/∑QₖXₖ
Ψₖₘ = exp(-(Aₖₘ + BₖₘT + CₖₘT²)/T)
```
## References
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. [doi:10.1021/i260064a004](https://doi.org/10.1021/i260064a004)
2. Weidlich, U.; Gmehling, J. A modified UNIFAC model. 1. Prediction of VLE, hE, and.gamma..infin. Ind. Eng. Chem. Res. 1987, 26, 1372–1381.
"""
UNIFAC

function UNIFAC(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

    groups = GroupParam(components, ["Activity/UNIFAC/UNIFAC_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)

    params = getparams(groups, ["Activity/UNIFAC/UNIFAC_like.csv", "Activity/UNIFAC/UNIFAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    R  = params["R"]
    Q  = params["Q"]
    _puremodel = init_puremodel(puremodel,groups.components,pure_userlocations,verbose)
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String["10.1021/i260064a004"]
    cache = UNIFACCache(groups,packagedparams)
    model = UNIFAC(groups.components,groups,packagedparams,_puremodel,references,cache)
    return model
end

function recombine_impl!(model::UNIFACModel)
    recombine_unifac_cache!(model.unifac_cache,model.groups,model.params)
    recombine!(model.puremodel)
    return model
end

function activity_coefficient(model::UNIFACModel,V,T,z)
    return exp.(@f(lnγ_comb)+ @f(lnγ_res))
end

function lnγ_comb(model::UNIFACModel,V,T,z)
    x = z ./ sum(z)
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    Φ = r/dot(x,r)
    Φ_p = q_p/dot(x,q_p)
    θ = q/dot(x,q)
    lnγ_comb = @. log(Φ_p)+(1-Φ_p)-5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end

function lnγ_SG(model::UNIFACModel,V,T,z)

    x = z ./ sum(z)

    r =model.unifac_cache.r
    q =model.unifac_cache.q

    Φ = r/dot(x,r)
    θ = q/dot(x,q)
    lnγ_SG = @. -5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_SG
end

function lnγ_res(model::UNIFACModel,V,T,z)
    v  = model.groups.n_flattenedgroups
    _ψ = @f(Ψ)
    lnΓ_ = @f(lnΓ,_ψ)
    lnΓi_ = @f(lnΓi,_ψ)
    lnγ_res_ =  [sum(v[i][k].*(lnΓ_[k].-lnΓi_[i][k]) for k ∈ @groups) for i ∈ @comps]
    return lnγ_res_
end

function lnΓ(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    x = z ./ sum(z)
    X = sum(v[i][:]*x[i] for i ∈ @comps) ./ sum(sum(v[i][k]*x[i] for k ∈ @groups) for i ∈ @comps)
    θ = X.*Q / dot(X,Q)
    lnΓ_ = Q.*(1 .-log.(sum(θ[m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[m]*ψ[:,m]./sum(θ[n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups))
    return lnΓ_
end

function lnΓi(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    ψ = @f(Ψ)
    X = [v[i][:] ./ sum(v[i][k] for k ∈ @groups) for i ∈ @comps]
    θ = [X[i][:].*Q ./ sum(X[i][n]*Q[n] for n ∈ @groups) for i ∈ @comps]
    lnΓi_ = [Q.*(1 .-log.(sum(θ[i][m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[i][m]*ψ[:,m]./sum(θ[i][n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups)) for i ∈ @comps]
    return lnΓi_
end

function Ψ(model::UNIFACModel,V,T,z)
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    return @. exp(-(A+B*T+C*T^2)/T)
end

function excess_g_SG(model::UNIFACModel,p,T,z)
    lnγ = lnγ_SG(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function excess_g_res(model::UNIFACModel,p,T,z)
    lnγ = lnγ_res(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end