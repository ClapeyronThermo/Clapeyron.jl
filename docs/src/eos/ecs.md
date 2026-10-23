```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["ecs.md"]
```

## Potential models

Potential equation of state models are a particular subtype of models that represent a 

## LJ Potentials

```@docs
Clapeyron.Thol2016LJ
Clapeyron.Johnson1993LJ
Clapeyron.KolafaNezbedaLJ
```

## LJ-TS Potentials

```@docs
Clapeyron.PeTSLJ
```

## Mie Potentials

```@docs
Clapeyron.FeANNMie
Clapeyron.Pohl2023Mie6
```

## Corresponding states models

Corresponding states models inheriting the `CorrespondingStatesModel` in Clapeyron use the following API:

```julia
scaling = Clapeyron.scaling_model(model::CorrespondingStatesModel) #get inner scaling model
reference = Clapeyron.reference_model(model::CorrespondingStatesModel) #get inner reference model
f,h = Clapeyron.shape_factors(scaling::MyScalingModel,reference,T,z) #get temperature and volume scaling factors, as a method of the scaling model.
```

Clapeyron provides specialized methods that reroute the evaluation of some properties to the reference model first, and then perform the scaling, instead of creating shape factors multiple times.

## Corresponding states API

```@docs
Clapeyron.CorrespondingStatesModel
Clapeyron.scaling_model
Clapeyron.reference_model
Clapeyron.shape_factors
```

## Scaling models

```@docs
Clapeyron.NoECSScaling
Clapeyron.SPUNGScaling
Clapeyron.LJScaling
Clapeyron.ChainLJScaling
```

## Corresponding States models

```@docs
Clapeyron.ECS
Clapeyron.SPUNG
Clapeyron.KolafaNezbeda
Clapeyron.PeTS
Clapeyron.LJRef
```
