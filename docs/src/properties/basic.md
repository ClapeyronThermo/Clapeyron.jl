```@meta
CurrentModule = Clapeyron
```

## Contents

```@contents
Pages = ["basic.md"]
```

## Index

```@index
Pages = ["basic.md"]
```

## Primitive functions

Almost all models in Clapeyron based on helmholtz free energy have at least one of the following functions defined:

```@docs
Clapeyron.eos
Clapeyron.eos_res
Clapeyron.idealmodel
Clapeyron.a_res
```

## Automatic Differenciation functions

All bulk properties in `Clapeyron` are calculated via a combination of these Automatic Differenciation Primitives over [`eos`](@ref) or [`eos_res`](@ref)
 
```@docs
Clapeyron.∂f∂T
Clapeyron.∂f∂V
Clapeyron.∂f
Clapeyron.p∂p∂V
Clapeyron.∂2f
Clapeyron.∂2p
Clapeyron.f_hess
Clapeyron.∂²³f
```

## Thermodynamic Method Dispatch types
```@docs
Clapeyron.ThermodynamicMethod
Clapeyron.SaturationMethod
Clapeyron.BubblePointMethod
Clapeyron.DewPointMethod
Clapeyron.TPFlashMethod
```