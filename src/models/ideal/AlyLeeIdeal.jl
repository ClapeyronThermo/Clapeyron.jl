struct AlyLeeIdealParam <: EoSParam
    Mw::SingleParam{Float64}
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    E::SingleParam{Float64}
    F::SingleParam{Float64}
    G::SingleParam{Float64}
    H::SingleParam{Float64}
    I::SingleParam{Float64}
end

abstract type AlyLeeIdealModel <: IdealModel end
@newmodelsimple AlyLeeIdeal AlyLeeIdealModel AlyLeeIdealParam

"""
    AlyLeeIdeal <: AlyLeeIdealModel
    AlyLeeIdeal(components::Array{String,1}; 
    userlocations::Array{String,1}=String[], 
    verbose=false)

## Input parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
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

## References

1. Aly, F. A., & Lee, L. L. (1981). Self-consistent equations for calculating the ideal gas heat capacity, enthalpy, and entropy. Fluid Phase Equilibria, 6(3–4), 169–179. [doi:10.1016/0378-3812(81)85002-9](https://doi.org/10.1016/0378-3812(81)85002-9)
"""

function AlyLeeIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose, ignore_missing_singleparams = ["B","C","D","E","F","G","H","I"])
    Mw = params["Mw"]
    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]
    E = params["E"]
    F = params["F"]
    G = params["G"]
    H = params["H"]
    I = params["I"]
    packagedparams = AlyLeeIdealParam(Mw,A,B,C,D,E,F,G,H,I)
    references = ["10.1016/0378-3812(81)85002-9"]
    #
    return AlyLeeIdeal(packagedparams;references,verbose)
end


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
    Tinv = one(T)/T

    @inbounds for i ∈ @comps
        #we suppose ρc = Tc = 1
        δi = ρ
        τi = Tinv
        Ai,Bi,Ci,Di,Ei,Fi,Gi,Hi,Ii = A[i],B[i],C[i],D[i],E[i],F[i],G[i],H[i],I[i]
        #integrate constant:
        #Tc,T0 = 1.0,298.15
        #τ0 = T0/Tc
        a₁ = Ai*-4.697596715569115 #(1 - log(τ0))
        a₂ = -Ai*298.15 #T0/Tc
        c₀ = (Ai - 1)
        ai = a₁ + a₂*τi + c₀*log(τi)
        iszero(Bi) || (ai += Bi*log(abs(sinh(Ci*τi))))
        iszero(Di) || (ai -= Di*log(cosh(Ei*τi)))
        iszero(Fi) || (ai += Fi*log(abs(sinh(Gi*τi)))) 
        iszero(Hi) || (ai -= Gi*log(cosh(Ii*τi)))
        res += z[i]*(ai + log(δi) + log(z[i]) - lnΣz)
    end
    return res/Σz
end

#=
ci = polycoeff[i]
n = length(ci)
c0 = first(ci)
cii = last(ci,n-1)
div1 = NTuple{n,Int}(1:n)
div2 = NTuple{n-1,Int}(1:n)
R̄⁻¹= 1/R̄
pol1 = ci ./ div1
pol2 = cii ./ div2
lnV = log(z[i]/V)
lnT = (1 - c0*R̄⁻¹)*(log(T) - log(298))

H = (c0*T - 298*c0*298/T)*R̄⁻¹
res += z[i]*(lnV+lnT+H)


=#

export AlyLeeIdeal

