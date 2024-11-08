@newmodelsimple Soave2019Alpha SoaveAlphaModel SimpleAlphaParam
export Soave2019Alpha

"""
    Soave2019Alpha <: SoaveAlphaModel

    Soave2019Alpha(components::Vector{String};
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description
Cubic alpha `(α(T))` model. Updated m(ω) correlations for `PR` and `SRK` with better results for heavy molecules.
```
αᵢ = (1+mᵢ(1-√(Trᵢ)))^2
Trᵢ = T/Tcᵢ
mᵢ = 0.37464 + 1.54226ωᵢ - 0.26992ωᵢ^2
```
where, for Peng-robinson:
```
mᵢ = 0.3919 + 1.4996ωᵢ - 0.2721ωᵢ^2 + 0.1063ωᵢ^3
```
and, for Redlich-Kwong:
```
mᵢ = 0.4810 + 1.5963ωᵢ - 0.2963ωᵢ^2 + 0.1223ωᵢ^3
```

## Model Construction Examples
```
# Using the default database
alpha = Soave2019Alpha("water") #single input
alpha = Soave2019Alpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = Soave2019Alpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/acentric.csv"])

# Passing parameters directly
alpha = Soave2019Alpha(["neon","hydrogen"];userlocations = (;acentricfactor = [-0.03,-0.21]))
```

## References
1. Pina-Martinez, A., Privat, R., Jaubert, J.-N., & Peng, D.-Y. (2019). Updated versions of the generalized Soave α-function suitable for the Redlich-Kwong and Peng-Robinson equations of state. Fluid Phase Equilibria, 485, 264–269. [doi:10.1016/j.fluid.2018.12.007](https://doi.org/10.1016/j.fluid.2018.12.007)
"""
Soave2019Alpha

default_locations(::Type{Soave2019Alpha}) = critical_data()
default_references(::Type{Soave2019Alpha}) = ["10.1016/j.fluid.2018.12.007"]
@inline α_m(model::PRModel,::Soave2019Alpha) = (0.3919,1.4996,-0.2721,0.1063)
@inline α_m(model::RKModel,::Soave2019Alpha) = (0.4810,1.5963,-0.2963,0.1223)