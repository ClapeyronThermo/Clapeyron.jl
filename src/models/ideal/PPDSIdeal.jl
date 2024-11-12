struct PPDSIdealParam <: EoSParam
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    E::SingleParam{Float64}
    F::SingleParam{Float64}
    G::SingleParam{Float64}
    reference_state::ReferenceState
end

abstract type PPDSIdealModel <: IdealModel end
@newmodelsimple PPDSIdeal PPDSIdealModel PPDSIdealParam
export PPDSIdeal
"""
    PPDSIdeal <: PPDSIdealModel

    PPDSIdeal(components;
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

## Description

PPDS Ideal Model:

```
Cpᵢ(T)/R = B + (C - B)y²[1 + (y − 1)(D + Ey + Fy² + Gy³)]
y = T/(A + T)
```

## Model Construction Examples
```
# Using the default database
idealmodel = PPDSIdeal("water") #single input
idealmodel = PPDSIdeal(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = PPDSIdeal(["neon","hydrogen"]; userlocations = ["path/to/my/db","alylee.csv"])

# Passing parameters directly
idealmodel = PPDSIdeal(["water","carbon dioxide"];
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

1. Gmehling, J., Kleiber, M., Kolbe, B., & Rarey, J. (2019). Chemical thermodynamics for process simulation (2nd ed.). Berlin, Germany: Blackwell Verlag.
"""
PPDSIdeal
default_locations(::Type{PPDSIdeal}) = ["ideal/PPDSIdeal.csv"]
default_references(::Type{PPDSIdeal}) = ["978-3-527-34325-6"] #TODO: find original source of the equation.

function a_ideal_T(model::PPDSIdealModel,T,z)
    #we transform from AlyLee terms to GERG2008 terms.

    _A = model.params.A.values
    _B = model.params.B.values
    _C = model.params.C.values
    _D = model.params.D.values
    _E = model.params.E.values
    _F = model.params.F.values
    _G = model.params.G.values

    Σz = sum(z)
    res = zero(T+first(z))
    τi = one(T)/T
    Tr = one(T)
    logτi = log(τi)
    @inbounds for i ∈ @comps
        #we suppose ρc = Vc = 1
        zi = z[i]
        A,B,C,D,E,F,G = _A[i],_B[i],_C[i],_D[i],_E[i],_F[i],_G[i]


            #helmholtz formulation in 10.1007/s10765-024-03360-0, eq 23
            #we suppose CI == CII == 0, ther reference_state field takes care of that.
            η = (B - C)
            λ = G + F + E + D
            Ā = A*τi + Tr
            logĀ = log(Ā)

            ai = η*(λ + 2)*(Ā*logĀ - A*τi)/Tr
            #@show A
            ai += A*(C - B)*(λ + 2)*τi*(logτi - 1.0)/Tr
            #@show ai
            ai += -η*(λ + 1)*logĀ
            #@show ai
            ai += (C - 1)*logτi
            #@show ai
            y = Tr/Ā
            coeffs = (zero(λ), 0.5*λ, (λ - D)/6, (G + F)/12, G/20)
            ai += η*evalpoly(y,coeffs)
            res += zi*ai
    end
    return res/Σz
end
