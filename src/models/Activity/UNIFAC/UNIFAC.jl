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
    icomponents::UnitRange{Int}
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

    UNIFAC(components::Vector{String};
    puremodel = PR,
    userlocations = String[], 
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

The Combinatorial part corresponds to an GC-averaged modified UNIQUAC model. The residual part iterates over groups instead of components.

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
1. Fredenslund, A., Gmehling, J., Michelsen, M. L., Rasmussen, P., & Prausnitz, J. M. (1977). Computerized design of multicomponent distillation columns using the UNIFAC group contribution method for calculation of activity coefficients. Industrial & Engineering Chemistry Process Design and Development, 16(4), 450–462. doi:10.1021/i260064a004
2. Weidlich, U.; Gmehling, J. A modified UNIFAC model. 1. Prediction of VLE, hE, and.gamma..infin. Ind. Eng. Chem. Res. 1987, 26, 1372–1381.
"""
UNIFAC

function UNIFAC(components::Vector{String};
    puremodel = PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false)
    
    groups = GroupParam(components, ["Activity/UNIFAC/UNIFAC_groups.csv"]; verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/UNIFAC_like.csv", "Activity/UNIFAC/UNIFAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    R  = params["R"]
    Q  = params["Q"]
    icomponents = 1:length(components)
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String["10.1021/i260064a004"]
    cache = UNIFACCache(groups,packagedparams)
    model = UNIFAC(components,icomponents,groups,packagedparams,_puremodel,references,cache)
    return model
end

function excess_g_SG(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    θm =  dot(q,z)*invn
    G_sg = _0
    for i ∈ @comps
        xi = z[i]*invn
        Φi = r[i]/Φm
        θi = q[i]/θm
        G_sg += 5*q[i]*xi*log(θi/Φi)
    end
    return n*G_sg
end

function excess_g_comb(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    Φpm = dot(q_p,z)*invn
    θm =  dot(q,z)*invn
    G_comb = _0
    for i ∈ @comps
        xi = z[i]*invn
        Φi = r[i]/Φm
        Φpi = q_p[i]/Φpm
        θi = q[i]/θm
        G_comb += xi*log(Φpi) + 5*q[i]*xi*log(θi/Φi)
    end
    return n*G_comb
end

function excess_g_res(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(T+first(z))
    Q = model.params.Q.values
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    invT = 1/T
    mi = group_matrix(model.groups)
    m̄ = dot(z,model.unifac_cache.m)
    m̄inv = 1/m̄
    X = m̄inv*mi*z
    #X = [dot(z,mi_i)*m̄inv for mi_i in eachcol(mi)]
    θpm = dot(X,Q)
    G_res = _0
    for i ∈ @groups
        q_pi = Q[i]
        ∑θpτ = _0
        for j ∈ @groups
            θpj = Q[j]*X[j]/θpm
            τji = exp(-evalpoly(T,(A[j,i],B[j,i],C[j,i]))*invT)
            ∑θpτ += θpj*τji
        end
        G_res += q_pi*X[i]*log(∑θpτ)
    end
    return -m̄*G_res
end

function excess_gibbs_free_energy(model::UNIFACModel,p,T,z)
    g_comp = excess_g_comb(model,p,T,z)
    g_res = excess_g_res(model,p,T,z)
    return (g_comp+g_res)*R̄*T 
end
