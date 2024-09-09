struct AlyLeeIdealParam <: EoSParam
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    E::SingleParam{Float64}
    F::SingleParam{Float64}
    G::SingleParam{Float64}
    H::SingleParam{Float64}
    I::SingleParam{Float64}
    reference_state::ReferenceState
end

abstract type AlyLeeIdealModel <: IdealModel end
@newmodelsimple AlyLeeIdeal AlyLeeIdealModel AlyLeeIdealParam

"""
    AlyLeeIdeal <: AlyLeeIdealModel

    AlyLeeIdeal(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters

- `A`: Single Parameter (`Float64`) - Model Coefficient
- `B`: Single Parameter (`Float64`) - Model Coefficient
- `C`: Single Parameter (`Float64`) - Model Coefficient
- `D`: Single Parameter (`Float64`) - Model Coefficient
- `E`: Single Parameter (`Float64`) - Model Coefficient
- `F`: Single Parameter (`Float64`) - Model Coefficient
- `G`: Single Parameter (`Float64`) - Model Coefficient
- `H`: Single Parameter (`Float64`) - Model Coefficient
- `I`: Single Parameter (`Float64`) - Model Coefficient


## Description

Aly-Lee Ideal Model (extended):

```
Cpᵢ(T)/R = A + B(CT⁻¹/sinh(CT⁻¹))² + D(ET⁻¹/cosh(ET⁻¹))² + F(GT⁻¹/sinh(GT⁻¹))² + H(IT⁻¹/cosh(IT⁻¹))²
```

## Model Construction Examples
```
# Using the default database
idealmodel = AlyLeeIdeal("water") #single input
idealmodel = AlyLeeIdeal(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = AlyLeeIdeal(["neon","hydrogen"]; userlocations = ["path/to/my/db","alylee.csv"])

# Passing parameters directly
idealmodel = AlyLeeIdeal(["water","carbon dioxide"];
                        userlocations = (A = [4.004, 3.5],
                        B = [0.01, 2.044],
                        C = [268.8, 919.3],
                        D = [0.99, -1.06],
                        E = [1141.4, -865.1],
                        F = [3.07, 2.034],
                        G = [2507.37, 483.55],
                        H = [0.0, 0.0139],
                        I = [0.0, 341.11])
                        )
```

## References

1. Aly, F. A., & Lee, L. L. (1981). Self-consistent equations for calculating the ideal gas heat capacity, enthalpy, and entropy. Fluid Phase Equilibria, 6(3–4), 169–179. [doi:10.1016/0378-3812(81)85002-9](https://doi.org/10.1016/0378-3812(81)85002-9)
"""
AlyLeeIdeal
default_locations(::Type{AlyLeeIdeal}) = ["ideal/AlyLeeIdeal.csv"]
default_references(::Type{AlyLeeIdeal}) = ["10.1016/0378-3812(81)85002-9"]

function a_ideal(model::AlyLeeIdealModel,V,T,z=SA[1.0])
    #we transform from AlyLee terms to GERG2008 terms.

    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    D = model.params.D.values
    E = model.params.E.values
    F = model.params.F.values
    G = model.params.G.values
    H = model.params.H.values
    I = model.params.I.values
    
    Σz = sum(z)
    res = zero(V+T+first(z))
    ρ = Σz/V
    lnΣz = log(Σz)
    τi = one(T)/T
    logτi = log(τi)
    @inbounds for i ∈ @comps
        #we suppose ρc = Tc = 1
        δi = ρ
        Ai,Bi,Ci,Di,Ei,Fi,Gi,Hi,Ii = A[i],B[i],C[i],D[i],E[i],F[i],G[i],H[i],I[i]
        #integrate constant:
        #Tc,T0 = 1.0,298.15
        a₁,a₂,c₀ = _Cp0_constant_parse(Ai,1.0,298.15)
        #a₁ = Ai*-4.697596715569115 #(1 - log(τ0))
        #a₂ = -Ai*298.15 #T0/Tc
        c₀ -= 1
        ai = a₁ + a₂*τi + c₀*logτi
        zi = z[i]
        ni = (Bi,Di,Fi,Hi)
        vi = (Ci,Ei,Gi,Ii)
        ai += term_a0_gerg2008(τi,logτi,zero(τi),ni,vi)
        res += xlogx(zi)
        res += zi*(ai + log(δi) - lnΣz)
    end
    return res/Σz - 1
end

export AlyLeeIdeal