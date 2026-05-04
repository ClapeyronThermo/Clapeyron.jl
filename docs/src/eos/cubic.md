```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["cubic.md"]
```

# Cubic Models

All cubic models in `Clapeyron.jl` follow a common evaluation order:

```julia
function CubicModel(args...)
    # get params for database, initialize other models, etc.
    recombine!(model)  # we calculate the mixing rules, caches for the translation models if necessary, etc.
end

function cubic_ab(model::CubicModel,V,T,z=SA[1.0])
    invn2 = (one(n)/n)^2
    a = model.params.a.values
    b = model.params.b.values
    α = α_function(model,V,T,z,model.alpha)
    c = translation(model,V,T,z,model.translation)
    ā,b̄,c̄ = mixing_rule(model,V,T,z,model.mixing,α,a,b,c)
    return ā, b̄, c̄
end

function a_res(model::CubicModel,V,T,z,data = (sum(z),cubic_ab(model,V,T,z)))
    n, ā, b̄, c̄ = data
    # depends on the specific EoS
    return result
end
```

- A *Mixing Rule Model* creates `aᵢⱼ` and `bᵢⱼ` from the critical temperature, critical pressure and a matrix of pair coefficients.
- An *Alpha Model* creates a vector of `αᵢ(T)` values.
- A *Translation Model* creates a vector of `cᵢ` values.
- The same Mixing rule, given `aᵢⱼ`, `bᵢⱼ`, `αᵢ(T)` and `cᵢ` returns the the mixture values of `ā`, `b̄` and `c̄` that are then used by the corresponding cubic model.
  A Mixing Rule can contain activity models to participate in the mixing (for example, Huron–Vidal rules).

## Common Definitions

```@docs
Clapeyron.ab_premixing
Clapeyron.mixing_rule
```

## van der Walls Models

```@docs
Clapeyron.vdW
Clapeyron.Clausius
Clapeyron.Berthelot
```

## Redlich-Kwong Models

```@docs
Clapeyron.RK
Clapeyron.SRK
Clapeyron.PSRK
Clapeyron.tcRK
```

## Peng-Robinson Models

```@docs
Clapeyron.PR
Clapeyron.PR78
Clapeyron.cPR
Clapeyron.VTPR
Clapeyron.TVTPR
Clapeyron.UMRPR
Clapeyron.tcPR
Clapeyron.tcPRW
Clapeyron.EPPR78
Clapeyron.QCPR
```

## Patel-Teja Models

```@docs
Clapeyron.PatelTeja
Clapeyron.PatelTejaHayen
Clapeyron.PTV
Clapeyron.YFR
```

## Other cubic Models

```@docs
Clapeyron.RKPR
Clapeyron.KU
```

## Alpha `(α(T))` Models

### Alpha (α(T)) Models: Family, Default Usage, and Compatibility with Cubic EoS

| Alpha Model | Alpha Family | Default Model(s) | Compatible with PR? | Compatible with RK? | Compatible with Other Cubic EoS? |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`NoAlpha`** (α=1) | — | `vdW` | Yes | Yes | Yes |
| **`ClausiusAlpha`** | — | `Clausius`, `Berthelot` | Yes | Yes | Yes |
| **`RKAlpha`** | — | `RK` | Yes | Yes | Yes |
| **`SoaveAlpha`** | Soave | `SRK` | Yes | Yes | `vdW` |
| **`Soave2019Alpha`** | Soave | None | Yes | Yes | `vdW` |
| **`PRAlpha`** | Soave | `PR` | Yes | No | No |
| **`PR78Alpha`** | Soave | `PR78`, `EPPR78` | Yes | No | No |
| **`PatelTejaAlpha`** | Soave | `PatelTeja` | No | No | No |
| **`PTVAlpha`** | Soave | `PTV` | No | No | specific correlation for `PTV` |
| **`MTAlpha`** | Soave | `UMRPR` | Yes | No | Yes (`vdW`) |
| **`BMAlpha`** | Soave (with extrapolation when T > Tc) | None | Yes | Yes | Yes (estimated)² |
| **`LeiboviciAlpha`** | Soave | None | Yes² | Yes² | Yes² |
| **`TwuAlpha`** / `Twu91Alpha` | Twu | `VTPR`, `TVTPR` | Yes¹ | Yes¹ | Yes¹ |
| **`tcTwuAlpha`** | Twu | Provides estimation for PR/RK | Yes¹ | Yes¹ | Yes (with custom parameters) |
| **`RKPRAlpha`** | — | `RKPR` | Yes (estimated)² | Yes (estimated)² | Yes (estimated)², specific correlation for `RKPR`|
| **`KUAlpha`** | Soave-ish (with extrapolation when T > Tc) | `KU` | No³ | No³ | No³ |
| **`YFRAlpha`** | — | `YFR` | No³ | No³ | specific correlation for `YFR` |
| **`PTHAlpha`** | — | `PatelTejaHayen` (only) | No³ | No³ | specific correlation for `PatelTejaHayen` |
| **`MathiasCopemanAlpha`** | Mathias-Copeman | None | Only with custom params | Only with custom params | Only with custom params |
| **`MC3PRAlpha`** | Mathias-Copeman | Provides estimation for PR | Yes (estimated) | Only with custom params | Only with custom params |
| **`CPAAlpha`** / `sCPAAlpha` | Soave | CPA, sCPA | Only with custom params | Only with custom params | Only with custom params |

**Footnotes:**

1. **Twu Family**: Parameters (L, M, N) can be automatically estimated for PR and RK when using `tcTwuAlpha`.  For base `TwuAlpha` with other cubics, estimation may be available depending on the EoS (e.g., with `tcPR`, `tcRK`, `QCPR`, `VTPR`); otherwise custom parameters are required.

2. **`LeiboviciAlpha`** is a unified Soave-type model that generates the `m(ω)` parameter based on the specific cubic EoS constants `Δ₁` and `Δ₂`. Consequently, it can be used with **any** cubic equation of state without requiring EoS-specific parameter correlations. Because of its generality, this model is used to provide estimation schemes for some alpha functions.

3. **Model-Specific Alphas**: `KUAlpha`, `YFRAlpha`, and `PTHAlpha` are tightly coupled to their respective EoS (`KU`, `YFR`, `PatelTejaHayen`) with unique functional forms and/or parameter correlations. They are not designed or documented for use with PR, RK, or other cubic families.

```@docs
Clapeyron.α_function
Clapeyron.NoAlpha
Clapeyron.ClausiusAlpha
Clapeyron.RKAlpha
Clapeyron.SoaveAlpha
Clapeyron.Soave2019Alpha
Clapeyron.PRAlpha
Clapeyron.PR78Alpha
Clapeyron.CPAAlpha
Clapeyron.sCPAAlpha
Clapeyron.MTAlpha
Clapeyron.BMAlpha
Clapeyron.TwuAlpha
Clapeyron.Twu91Alpha
Clapeyron.Twu88Alpha
Clapeyron.tcTwuAlpha
Clapeyron.PTHAlpha
Clapeyron.PatelTejaAlpha
Clapeyron.PTVAlpha
Clapeyron.KUAlpha
Clapeyron.RKPRAlpha
Clapeyron.LeiboviciAlpha
Clapeyron.MathiasCopemanAlpha
Clapeyron.MC3PRAlpha
```

## Volume Translation Models

```@docs
Clapeyron.translation
Clapeyron.NoTranslation
Clapeyron.ConstantTranslation
Clapeyron.ClausiusTranslation
Clapeyron.PenelouxTranslation
Clapeyron.RackettTranslation
Clapeyron.TVTPRTranslation
Clapeyron.MTTranslation
Clapeyron.BaledTranslation
```

## Mixing Rule Models

```@docs
Clapeyron.vdW1fRule
Clapeyron.KayRule
Clapeyron.HVRule
Clapeyron.MHV1Rule
Clapeyron.MHV2Rule
Clapeyron.LCVMRule
Clapeyron.WSRule
Clapeyron.modWSRule
Clapeyron.VTPRRule
Clapeyron.PSRKRule
Clapeyron.UMRRule
Clapeyron.QCPRRule
Clapeyron.PPR78Rule
Clapeyron.gErRule
```
