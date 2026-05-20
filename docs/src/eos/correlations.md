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

## Saturation models and types

```@docs
Clapeyron.SaturationModel
Clapeyron.SaturationCorrelation
Clapeyron.LeeKeslerSat
Clapeyron.DIPPR101Sat
Clapeyron.AntoineEqSat
```

### ML-based saturation models

!!! note
    The following methods are provided by the companion package `MLThermoProperties`.

```@docs ;canonical = false
MLThermoProperties.GRAPPA
```

# Liquid Volume Correlations

Liquid volume correlations are any `EoSModel` that are subtypes of `LiquidVolumeModel`.

Most liquid volume correlations are only defined by a volume equation, but some models are defined in terms of the gibbs free energy.

```@docs
RackettLiquid
YamadaGunnLiquid
COSTALD
GrenkeElliottWater
HoltenWater
ZeroLiquid
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
Intended to be used in conjunction with a liquid model within a [CompositeModel](@ref).

```@docs
SolidHfus
SolidKs
IAPWS06
JagerSpanSolidCO2
```
