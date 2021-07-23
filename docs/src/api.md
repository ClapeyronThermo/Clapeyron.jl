

## Contents

```@contents
Pages = ["api.md"]
```

## Index

```@index
Pages = ["api.md"]
```

## Volume-Temperature Based Properties

```@docs
Clapeyron.second_virial_coefficient
Clapeyron.pip
```
## Pressure-Temperature Based Properties

```@docs
Clapeyron.volume
Clapeyron.compressibility_factor
```
## Single component properties

```@docs
Clapeyron.sat_pure
Clapeyron.enthalpy_vap
Clapeyron.crit_pure
Clapeyron.acentric_factor
```

## Automatic Differenciation functions
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

## Initial guess functions
```@docs
Clapeyron.lb_volume
Clapeyron.T_scale
Clapeyron.p_scale
Clapeyron.x0_volume
Clapeyron.x0_volume_liquid
Clapeyron.x0_volume_gas
Clapeyron.volume_virial
Clapeyron.x0_sat_pure
Clapeyron.x0_crit_pure
```


