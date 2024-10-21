```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["empiric.md"]
```

# Empiric Helmholtz Models

Empiric (or MultiParameter) models in Clapeyron are composed of three different, but interacting parts:

- Pure Fluid parameters
- Mixing volume and temperature
- Departure model

Pure Fluids are instantiated from CoolProp JSON files, via the [`SingleFluid`](@ref) struct.
In theory, any pure fluid should be supported.
Furthermore, there is support for using directly the fluids defined in the CoolProp library:

```julia
julia> SingleFluid("Ethanol")
ERROR: cannot found component file R113. Try loading the CoolProp library by loading it.
Stacktrace:
 ....
julia> using CoolProp #loads the CoolProp library and allows access to their JSON.
julia> SingleFluid("Ethanol")
MultiParameter Equation of state for Ethanol:
 Polynomial power terms: 6
 Exponential terms: 10
 Gaussian bell-shaped terms: 9
```

Multicomponent models are a collection of `SingleFluid` models + a mixing model + a departure model:

```julia
julia> model = GERG2008(["water","carbon dioxide"])
MultiFluid{EmpiricAncillary, AsymmetricMixing, EmpiricDeparture} with 2 components:
 "water"
 "carbon dioxide"
Contains parameters: Mw, Tc, Pc, Vc, acentricfactor, lb_volume

julia> model.pures
2-element Vector{SingleFluid{EmpiricAncillary}}:
 SingleFluid{EmpiricAncillary}("water")
 SingleFluid{EmpiricAncillary}("carbon dioxide")

julia> model.mixing
AsymmetricMixing with 2 components:
 "water"
 "carbon dioxide"
Contains parameters: gamma_T, gamma_v, beta_T, beta_v

julia> model.departure
EmpiricDeparture with 2 components:
 "water"
 "carbon dioxide"
Contains parameters: F, parameters
```

## Generic Models

```@docs
Clapeyron.SingleFluid
Clapeyron.SingleFluidIdeal
Clapeyron.MultiFluid
Clapeyron.EmpiricIdeal
```

## SingleFluid Models

```@docs
Clapeyron.XiangDeiters
Clapeyron.IAPWS95
Clapeyron.PropaneRef
Clapeyron.TholLJ
Clapeyron.Ammonia2023
```

## MultiComponent models

```@docs
Clapeyron.LJRef
Clapeyron.GERG2008
Clapeyron.EOS_LNG
Clapeyron.HelmAct
```

## Mixing models

```@docs
Clapeyron.LinearMixing
Clapeyron.AsymmetricMixing
Clapeyron.LorentzBerthelotMixing
```

## Departure models

```@docs
Clapeyron.EmpiricDeparture
Clapeyron.departure_functions
Clapeyron.create_departure
Clapeyron.GEDeparture
Clapeyron.QuadraticDeparture
```
