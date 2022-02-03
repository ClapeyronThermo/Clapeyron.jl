## Contents

```@contents
Pages = ["models.md"]
```

## Index

```@index
Pages = ["models.md"]
```
## Ideal Models

All `Clapeyron.jl` models can be separated between an ideal and a residual contribution. The ideal contribution can be obtained via integration of the ideal isobaric heat capacity:

``\frac{A_{\mathrm{ideal}}}{Nk_\mathrm{B}T} = \sum_{i=1}^{N_{\mathrm{Component}}} x_i\left[\ln{\frac{\rho_i}{\rho_0}}
+ \frac{1}{Nk_\mathrm{B}T} \int_{T_0}^T \!\!C_{p,i}^0 dT + \frac{H_{0,i}}{Nk_\mathrm{B}T}- \frac{1}{Nk_{B}}\!\!\int_{T_0}^T \frac{C_{p,i}^0}{T} dT -\ln{\frac{T}{T_0}}-\frac{S_{0,i}}{Nk_\mathrm{B}} - 1\right]``

```@docs
Clapeyron.idealmodel
Clapeyron.BasicIdeal
Clapeyron.ReidIdeal
Clapeyron.JobackIdeal
Clapeyron.MonomerIdeal
Clapeyron.WalkerIdeal
```

## Cubic Models

All cubic models in `Clapeyron.jl` follow a common evaluation order:
```julia

function CubicModel(args...)
    #get params for database, initialize other models, etc
    #we need Tc,pc and kij for the generation of the matrices aᵢⱼ and bᵢⱼ
    a,b = a,b = ab_premixing(CubicModel,mixing_model,Tc,pc,k)
    #rest of the initialization, return model
end

function cubic_ab(model::CubicModel,V,T,z=SA[1.0])
    invn2 = (one(n)/n)^2
    a = model.params.a.values
    b = model.params.b.values
    α = α_function(model,V,T,z,model.alpha)
    c = translation(model,V,T,z,model.translation)
    ā,b̄,c̄ = mixing_rule(model,V,T,z,model.mixing,α,a,b,c)
    return ā ,b̄, c̄
end

function a_res(model::CubicModel,V,T,z)
    ā ,b̄, c̄ = cubic_ab(model,V,T,z)
    #depends on the specific EoS
    return result
end
```
- A *Mixing Rule Model* creates `aᵢⱼ` and `bᵢⱼ` from the critical temperature, critical pressure and a matrix of pair coefficients.

- An *Alpha Model* creates a vector of `αᵢ(T)` values

- A *Translation Model* creates a vector of `cᵢ` values

- The same Mixing rule, given `aᵢⱼ`, `bᵢⱼ`, `αᵢ(T)` and `cᵢ` returns the the mixture values of `ā`, `b̄` and `c̄` that are then used by the corresponding cubic model. a Mixing Rule can contain activity models to participate in the mixing (for example, Huron-Vidal rules)

### Models
```@docs
Clapeyron.vdW
Clapeyron.RK
Clapeyron.PR
Clapeyron.SRK
Clapeyron.PSRK
Clapeyron.PR78
Clapeyron.UMRPR
Clapeyron.VTPR
```

### Alpha `(α(T))` Models

```@docs
Clapeyron.α_function
Clapeyron.NoAlpha
Clapeyron.RKAlpha
Clapeyron.SoaveAlpha
Clapeyron.PRAlpha
Clapeyron.PR78Alpha
Clapeyron.CPAAlpha
Clapeyron.sCPAAlpha
Clapeyron.MTAlpha
Clapeyron.BMAlpha
Clapeyron.TwuAlpha
```

### Volume Translation Models

```@docs
Clapeyron.translation
Clapeyron.NoTranslation
Clapeyron.RackettTranslation
Clapeyron.PenelouxTranslation
Clapeyron.MTTranslation
```

### Mixing Rule Models
```@docs
Clapeyron.ab_premixing
Clapeyron.mixing_rule
Clapeyron.vdW1fRule
Clapeyron.KayRule
Clapeyron.HVRule
Clapeyron.MHV1Rule
Clapeyron.MHV2Rule
Clapeyron.LCVMRule
Clapeyron.WSRule
Clapeyron.VTPRRule
Clapeyron.PSRKRule
```
