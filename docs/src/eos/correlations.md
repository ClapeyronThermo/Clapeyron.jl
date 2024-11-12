```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["ideal.md"]
```

# Correlation Models

Correlation models are, as their name says, fitted equations that express one property of a compound.
They meant to be used in conjunction with other models (like Activity models that require a saturated liquid volume), or via a `CompositeModel`.
Because they only overload one property, the way to define a correlation is different than normal `EoSModel`s.

# Saturation Correlations

Saturation Correlations are any `EoSModel` that are subtypes of [`SaturationModel`](@ref) return `psat(T)` and the upper limit `(Tc,Pc)` pair.
To define saturation correlations, you need to overload:

```julia
function crit_pure(model::MySaturationModel <: SaturationModel)
    ...
    return (Tc,Pc,NaN)
end

function Clapeyron.saturation_pressure_impl(model::MySaturationModel <: SaturationModel,T,::SaturationCorrelation)
    ...
    return (psat,NaN,NaN)
end
```

## Saturation Models and Types

```@docs
Clapeyron.SaturationModel
Clapeyron.SaturationCorrelation
Clapeyron.LeeKeslerSat
Clapeyron.DIPPR101Sat
```

# Liquid Volume Correlations
Liquid Volume Correlations are any `EoSModel` that are subtypes of `LiquidVolumeModel`.
They return `volume(model,p,T,z, phase = :liquid)`.

```@docs
RackettLiquid
YamadaGunnLiquid
COSTALD
```

# Virial Models

Virial models are defined in terms of the second virial coefficient, `B(T,z)`.
The reduced residual Helmholtz energy is defined as:

``\frac{A_\mathrm{res}}{Nk_\mathrm{B}T} = \frac{B}{V}``,

To implement a virial model, it is necessary to overload `Clapeyron.second_virial_coefficient_impl(model::<:SecondVirialModel,T,z)`.

```@docs
AbbottVirial
TsonopoulosVirial
EoSVirial2
```

# Solid Models

Solid models provide simple approximations to the excess chemical potential in the solid phase.
Intended to be used in conjunction with a liquid model within a CompositeModel.

```@docs
SolidHfus
SolidKs
```
