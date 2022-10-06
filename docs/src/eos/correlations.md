```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["ideal.md"]
```
# Correlation Models

Correlation models are, as their name says, fitted equations that express one property of a compound. they meant to be used in conjunction with other models (like Activity models that require a saturated liquid volume), or via a `CompositeModel`. because they only overload one property, the way to define a correlation is different than normal `EoSModel`s

# Saturation Correlations

Saturation Correlations are any [`EoSModel`](@ref) that are subtypes of [`SaturationModel`](@ref). return `psat(T)` and the upper limit `(Tc,Pc)` pair. To define saturation correlations, you need to overload:

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
Liquid Volume Correlations are any [`EoSModel`](@ref) that are subtypes of [`LiquidVolumeModel`](@ref). 
They return `volume(model,p,T,z, phase = :liquid)`.

# Virial Models

Virial models are defined in terms of the second virial coefficient, `B(T,z)`. this allows for a fast calculation of `volume(model,p,T,z, phase = :gas)`. they cannot represent the liquid phase.

