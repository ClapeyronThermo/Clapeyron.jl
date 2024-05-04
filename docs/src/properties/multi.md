```@meta
CurrentModule = Clapeyron
```

## Contents

```@contents
Pages = ["multi.md"]
```

## Index

```@index
Pages = ["multi.md"]
```

## Multi component properties

```@docs
Clapeyron.bubble_pressure
Clapeyron.bubble_temperature
Clapeyron.dew_pressure
Clapeyron.dew_temperature
Clapeyron.azeotrope_pressure
Clapeyron.azeotrope_temperature
Clapeyron.LLE_pressure
Clapeyron.LLE_temperature
Clapeyron.VLLE_pressure
Clapeyron.VLLE_temperature
Clapeyron.crit_mix
Clapeyron.UCEP_mix
Clapeyron.UCST_mix
Clapeyron.gibbs_solvation
Clapeyron.cross_second_virial
Clapeyron.equivol_cross_second_virial
Clapeyron.sle_solubility
Clapeyron.slle_solubility
Clapeyron.eutectic_point
```

### Bubble/Dew methods
```@docs
Clapeyron.ChemPotBubblePressure
Clapeyron.FugBubblePressure
Clapeyron.ActivityBubblePressure
Clapeyron.ChemPotBubbleTemperature
Clapeyron.FugBubbleTemperature
Clapeyron.ActivityBubbleTemperature
Clapeyron.ChemPotDewPressure
Clapeyron.FugDewPressure
Clapeyron.ActivityDewPressure
Clapeyron.ChemPotDewTemperature
Clapeyron.FugDewTemperature
Clapeyron.ActivityDewTemperature
```

## Consistency and Stability

```@docs
Clapeyron.gibbs_duhem
Clapeyron.isstable
Clapeyron.VT_mechanical_stability
Clapeyron.VT_diffusive_stability
Clapeyron.VT_chemical_stability
Clapeyron.tpd
```

## TP Flash

```@docs
Clapeyron.tp_flash
Clapeyron.DETPFlash
Clapeyron.RRTPFlash
Clapeyron.MichelsenTPFlash
Clapeyron.MultiPhaseTPFlash
Clapeyron.MCFlashJL
Clapeyron.numphases
Clapeyron.supports_reduction
```
