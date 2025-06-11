abstract type LeiboviciAlphaModel <: GeneralizedSuaveAlphaModel end

const LeiboviciAlphaParam = SimpleAlphaParam

@newmodelsimple LeiboviciAlpha LeiboviciAlphaModel SimpleAlphaParam
export LeiboviciAlpha

"""
    LeiboviciAlpha <: LeiboviciAlphaModel

    LeiboviciAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters


- `acentricfactor`: Single Parameter (`Float64`)

## Description

Leibovici Cubic alpha `(α(T))` model.
Generalized soave model that works for all common cubic models.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = ∑aᵢωᵢ^(i-1)
aᵢ = ∑bᵢu₀^(i-1)
u₀ = (u + 2)* √(2/(1 + u + w)) - 2
u = - Δ1 - Δ2
w = Δ1*Δ2
```

## Model Construction Examples

```
# Using the default database
alpha = LeiboviciAlpha("water") #single input
alpha = LeiboviciAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = LeiboviciAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = LeiboviciAlpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

## References

1. Leibovici, C. F. (1994). A unified m(ω) relation for cubic equations of state. Fluid Phase Equilibria, 101, 1–2. [doi:10.1016/0378-3812(94)02603-3](https://doi.org/10.1016/0378-3812(94)02603-3)

"""
LeiboviciAlpha

default_locations(::Type{LeiboviciAlpha}) = critical_data()
default_references(::Type{LeiboviciAlpha}) = ["10.1016/0378-3812(94)02603-3"]

@inline function α_m(model::DeltaCubicModel,alpha_model::LeiboviciAlphaModel,i)
    coeff = α_m_leibovici(model)
    ω = alpha_model.params.acentricfactor.values[i]
    return evalpoly(ω,coeff)
end

function α_m_leibovici(model::DeltaCubicModel)
    Δ1,Δ2 = cubic_Δ(model,SA[1.0])
    u = - Δ1 - Δ2
    w = Δ1*Δ2
    u0 = (u + 2)*sqrt(2/(1 + u + w)) - 2
    α_m_leibovici(u0)
end

function α_m_leibovici(model::DeltaCubicModel,i)
    z = FillArrays.OneElement(i, length(model))
    Δ1,Δ2 = cubic_Δ(model,z)
    u = - Δ1 - Δ2
    w = Δ1*Δ2
    u0 = (u + 2)*sqrt(2/(1 + u + w)) - 2
    α_m_leibovici(u0)
end

function α_m_leibovici(u::Number)
    _0 = zero(u)
    A = Leibovici_consts
    a0 = evalpoly(u,A[1])
    a1 = evalpoly(u,A[2])
    a2 = evalpoly(u,A[3])
    a3 = evalpoly(u,A[4])
    a4 = evalpoly(u,A[5])
    return (a0,a1,a2,a3,a4)
end

const Leibovici_consts = ((0.61090441, -0.14983477, 0.0209346, -0.00227238, 0.00013033),
(1.68811978, -0.12822686, 0.01842279, -0.00218662, 0.00012884),
(-0.22275603, 0.03341062, -0.00574116, 0.00084555, -5.925e-5),
(0.03702243, -0.01018267, 0.00205506, -0.00035826, 2.773e-5),
(-0.00351123, 0.00158691, -0.00035556, 6.834e-5, -5.54e-6))
