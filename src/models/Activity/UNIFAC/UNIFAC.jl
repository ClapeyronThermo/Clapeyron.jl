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

default_locations(::Type{UNIFAC}) = ["Activity/UNIFAC/UNIFAC_like.csv", "Activity/UNIFAC/UNIFAC_unlike.csv"]
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

## List of groups available
|Name      |SMARTS                                           |Description            |
|----------|-------------------------------------------------|-----------------------|
|CH3       |[CX4;H3;!R]                                      |Methyl group           |
|CH2       |[CX4;H2;!R]                                      |Methylene group        |
|CH        |[CX4;H1;!R]                                      |                       |
|C         |[CX4;H0;!R]                                      |                       |
|CH2=CH    |[CX3;H2]=[CX3;H1]                                |                       |
|CH=CH     |[CX3;H1]=[CX3;H1]                                |                       |
|CH2=C     |[CX3;H2]=[CX3;H0]                                |                       |
|CH=C      |[CX3;H1]=[CX3;H0]                                |                       |
|ACH       |[cX3;H1]                                         |Aromatic CH            |
|AC        |[cX3;H0]                                         |                       |
|ACCH3     |[cX3;H0][CX4;H3]                                 |                       |
|ACCH2     |[cX3;H0][CX4;H2]                                 |                       |
|ACCH      |[cX3;H0][CX4;H1]                                 |                       |
|OH (P)    |[OH1;$([OH1][CX4H2])]                            |Primary alcohol        |
|CH3OH     |[CX4;H3][OX2;H1]                                 |Methanol               |
|H2O       |[OH2]                                            |Water                  |
|ACOH      |[cX3;H0;R][OX2;H1]                               |                       |
|CH3CO     |[CX4;H3][CX3](=O)                                |Methyl ketone          |
|CH2CO     |[CX4;H2][CX3](=O)                                |Methylene ketone       |
|CHO       |[CX3H1](=O)                                      |                       |
|CH3COO    |[CH3][CX3;H0](=[O])[O]                           |Acetate group          |
|CH2COO    |[CX4;H2][CX3](=[OX1])[OX2]                       |                       |
|HCOO      |[CX3;H1](=[OX1])[OX2]                            |                       |
|CH3O      |[CH3;!R][OH0;!R]                                 |Methoxy group          |
|CH2O      |[CH2;!R][OH0;!R]                                 |                       |
|CHO       |[C;H1;!R][OH0;!R]                                |                       |
|THF       |[CX4;H2;R][OX2;R][CX4;H2;R]                      |Tetrahydrofuran        |
|CH3NH2    |[CX4;H3][NX3;H2]                                 |Methylamine            |
|CH2NH2    |[CX4;H2][NX3;H2]                                 |                       |
|CHNH2     |[CX4;H1][NX3;H2]                                 |                       |
|CH3NH     |[CX4;H3][NX3;H1]                                 |                       |
|CH2NH     |[CX4;H2][NX3;H1]                                 |                       |
|CHNH      |[CX4;H1][NX3;H1]                                 |                       |
|CH3N      |[CX4;H3][NX3;H0]                                 |                       |
|CH2N      |[CX4;H2][NX3;H0]                                 |                       |
|ACNH2     |[c][NX3;H2]                                      |                       |
|AC2H2N    |[cX3H1][n][cX3H1]                                |                       |
|AC2HN     |[cX3H0][n][cX3H1]                                |                       |
|AC2N      |[cX3H0][n][cX3H0]                                |                       |
|CH3CN     |[CX4;H3][CX2]#[NX1]                              |Acetonitrile           |
|CH2CN     |[CX4;H2][CX2]#[NX1]                              |                       |
|COO       |[CX3,cX3](=[OX1])[OX2H0,oX2H0]                   |Ester group            |
|COOH      |[CX3](=[OX1])[O;H1]                              |Carboxylate group      |
|HCOOH     |[CX3;H1](=[OX1])[OX2;H1]                         |Formic acid            |
|CH2CL     |[CX4;H2;!$(C(Cl)(Cl))](Cl)                       |                       |
|CHCL      |[CX4;H1;!$(C(Cl)(Cl))](Cl)                       |                       |
|CCL       |[CX4;H0](Cl)                                     |                       |
|CH2CL2    |[CX4;H2;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)               |Dichloromethane        |
|CHCL2     |[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)               |                       |
|CCL2      |[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)               |                       |
|CHCL3     |[CX4;H1;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)|Chloroform             |
|CCL3      |[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)|                       |
|CCL4      |[CX4;H0]([Cl])([Cl])([Cl])([Cl])                 |Carbon tetrachloride   |
|ACCL      |[c][Cl]                                          |                       |
|CH3NO2    |[CX4;H3][NX3](=[OX1])([OX1])                     |Nitromethane           |
|CH2NO2    |[CX4;H2][NX3](=[OX1])([OX1])                     |                       |
|CHNO2     |[CX4;H1][NX3](=[OX1])([OX1])                     |                       |
|ACNO2     |[cX3][NX3](=[OX1])([OX1])                        |                       |
|CS2       |C(=S)=S                                          |Carbon disulfide       |
|CH3SH     |[SX2H][CX4;H3]                                   |Methanethiol           |
|CH2SH     |[SX2H][CX4;H2]                                   |                       |
|FURFURAL  |c1cc(oc1)C=O                                     |Furfural               |
|DOH       |[OX2;H1][CX4;H2][CX4;H2][OX2;H1]                 |                       |
|I         |[I]                                              |Iodine group           |
|BR        |[Br]                                             |Bromine group          |
|CH=-C     |[CX2;H1]#[CX2;H0]                                |                       |
|C=-C      |[CX2;H0]#[CX2;H0]                                |                       |
|DMSO      |[SX3H0](=[OX1])([CX4;H3])[CX4;H3]                |Dimethyl sulfoxide     |
|ACRY      |[CX3;H2]=[CX3;H1][CX2;H0]#[NX1;H0]               |Acrylate               |
|CL-(C=C)  |[$([Cl;H0]([C]=[C]))]                            |                       |
|C=C       |[CX3;H0]=[CX3;H0]                                |                       |
|ACF       |[cX3][F]                                         |                       |
|DMF       |[CX4;H3][N]([CX4;H3])[CX3;H1]=[O]                |Dimethylformamide      |
|HCON(CH2)2|[NX3]([CX4;H2])([CX4;H2])[CX3;H1](=[OX1])        |                       |
|CF3       |C(F)(F)F                                         |                       |
|CF2       |C(F)F                                            |                       |
|CF        |C(F)                                             |                       |
|CY-CH2    |[CH2;R]                                          |Cycloalkane group      |
|CY-CH     |[CH1;R]                                          |                       |
|CY-C      |[CH0;R]                                          |                       |
|OH (S)    |[OH1;$([OH1][CX4H1])]                            |Second hydroxyl group  |
|OH (T)    |[OH1;$([OH1][CX4H0])]                            |Tertiary hydroxyl group|
|CY-CH2O   |[CX4H2;R][OX2;R]                                 |                       |
|TRIOXAN   |[CX4H2;R][OX2;R]                                 |Trioxane               |
|CNH2      |[CX4H0][NH2]                                     |                       |
|NMP       |[OX1H0]=[C;R][NX3H0;R][CH3]                      |N-Methylpyrrolidone    |
|NEP       |[OX1H0]=[CH0X3;R][H0;R][CH2]                     |                       |
|NIPP      |[OX1H0;!R]=[CX3H0;R][NX3H0;R][C;!R]              |                       |
|NTBP      |[OX1H0;!R]=[CH0X3;R][NX3H0;R][CH0;!R]            |                       |
|CONH2     |[CX3H0](=[OX1H0])[NX3H2]                         |                       |
|CONHCH3   |[OX1H0;!R]=[CX3H0;!R][NH1X3;!R][CH3;!R]          |                       |
|CONHCH2)  |[CH2X4;!R][NH1X3;!R][CX3H0;!R]=[OX1H0;!R]        |                       |
"""
UNIFAC

function UNIFAC(components;
    puremodel = PR,
    userlocations = String[],
    group_userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

    groups = GroupParam(components, ["Activity/UNIFAC/UNIFAC_groups.csv"]; group_userlocations = group_userlocations, verbose = verbose)

    params = getparams(groups, default_locations(UNIFAC);
                        userlocations=userlocations,
                        asymmetricparams=["A","B","C"],
                        ignore_missing_singleparams=["A","B","C"],
                        verbose=verbose)
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
    if hasfield(typeof(model),:unifac_cache)
        recombine_unifac_cache!(model.unifac_cache,model.groups,model.params)
    end
    recombine!(model.puremodel)
    return model
end

function Ψ(model::UNIFACModel,V,T,z)
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    return @. exp(-(A+B*T+C*T^2)/T)
end

excess_g_comb(model::UNIFACModel,p,T,z=SA[1.0]) = excess_g_comb_dormund(model,p,T,z)

#flory-huggins(FH) + Staverman-Guggenheim (SG) contributions
function excess_g_comb_dormund(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    θm = dot(q,z)*invn
    Φpm = dot(q_p,z)*invn
    G_comb = _0
    for i ∈ @comps
        Φi = r[i]/Φm #technically xi[i]r[i]/Φm, but it gets cancelled out (log(θi/Φi))
        θi = q[i]/θm #technically xi[i]q[i]/θm, but it gets cancelled out (log(θi/Φi))
        Φpi = q_p[i]/Φpm #technically xi[i]q_p[i]/θpm, but it gets cancelled out (log(Φpi/xi))
        zi = z[i]
        G_comb += zi*log(Φpi) + 5*q[i]*zi*log(θi/Φi)
    end
    return G_comb*R̄*T
end

function excess_g_comb_original(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    θm = dot(q,z)*invn
    G_comb = _0
    for i ∈ @comps
        Φi = r[i]/Φm #technically xi[i]r[i]/Φm, but it gets cancelled out (log(θi/Φi))
        θi = q[i]/θm #technically xi[i]q[i]/θm, but it gets cancelled out (log(θi/Φi))
        zi = z[i]
        G_comb += zi*log(Φi) + 5*q[i]*zi*log(θi/Φi)
    end
    return G_comb*R̄*T
end


#just Staverman-Guggenheim (SG) contribution, used for UMR
function excess_g_SG(model::UNIFACModel,p,T,z)
    _0 = zero(eltype(z))
    r = model.unifac_cache.r
    q = model.unifac_cache.q
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    θm = dot(q,z)*invn
    G_SG = _0
    for i ∈ @comps
        Φi = r[i]/Φm #technically xi[i]r[i]/Φm, but it gets cancelled out (log(θi/Φi))
        θi = q[i]/θm #technically xi[i]q[i]/θm, but it gets cancelled out (log(θi/Φi))
        G_SG += 5*q[i]*z[i]*log(θi/Φi)
    end
    return G_SG*R̄*T
end

# https://github.com/thermotools/thermopack/blob/main/doc/memo/UNIFAC/unifac.pdf
function excess_g_res(model::UNIFACModel,p,T,z)
    _0 = zero(T + first(z))
    G_res = _0
    V = _0
    Ẽ = @f(Ψ)
    v = model.groups.n_flattenedgroups
    Q = model.params.Q.values
    ∑vikQk = [dot(Q,vi) for vi in v]
    #calculate Θ with the least amount of allocs possible
    X = group_matrix(model.groups)*z
    ∑XQ⁻¹ = 1/dot(X,Q)
    X .*= Q
    X .*= ∑XQ⁻¹
    Θ = X
    for i in @comps
        ∑QkΔΛk = _0
        vi = v[i]
        ∑vikQk⁻¹ = 1/∑vikQk[i] # 1/(sum(vik * Qk for k in groups))
        zi = z[i]
        #iszero(zi) && continue #causes problems with AD
        for k in @groups
            Λk = _0
            Λki = _0
            for j in @groups
                Ẽjk = Ẽ[j,k]
                Λk += Ẽjk*Θ[j]
                Θij = Q[j]*vi[j]
                Λki += Ẽjk*Θij
            end
            Λki = Λki*∑vikQk⁻¹
            Λk = log(Λk)
            Λki = log(Λki)
            ∑QkΔΛk += vi[k]*Q[k]*(Λk - Λki)
        end
        G_res +=zi*∑QkΔΛk
    end
    return -G_res*R̄*T
end

function excess_gibbs_free_energy(model::UNIFACModel,p,T,z)
    g_comb = excess_g_comb(model,p,T,z)
    g_res = excess_g_res(model,p,T,z)
    return g_comb+g_res
end

#=
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

function excess_g_comb(model::UNIFACModel,p,T,z)
    lnγ = lnγ_comb(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
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

function excess_g_res(model::UNIFACModel,p,T,z)
    lnγ = lnγ_res(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ 1:length(model))
end

function excess_g_SG(model::UNIFACModel,p,T,z)
    lnγ = lnγ_SG(model,p,T,z)
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end
=#
