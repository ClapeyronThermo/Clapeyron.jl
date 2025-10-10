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

## Bubble/Dew Points

```@docs
Clapeyron.bubble_pressure
Clapeyron.bubble_temperature
Clapeyron.dew_pressure
Clapeyron.dew_temperature
Clapeyron.gibbs_solvation
```

## Azeotropes, LLE and VLLE equilibria

```@docs
Clapeyron.azeotrope_pressure
Clapeyron.azeotrope_temperature
Clapeyron.LLE_pressure
Clapeyron.LLE_temperature
Clapeyron.LLE
Clapeyron.VLLE_pressure
Clapeyron.VLLE_temperature
```

## Critical points, spinodals and stability limits

```@docs
Clapeyron.crit_mix
Clapeyron.mechanical_critical_point
Clapeyron.spinodal_pressure
Clapeyron.spinodal_temperature
Clapeyron.edge_pressure
Clapeyron.edge_temperature
Clapeyron.UCEP_mix
Clapeyron.UCST_pressure
Clapeyron.UCST_temperature
```

## SLE Equilibria

```@docs
Clapeyron.sle_solubility
Clapeyron.sle_solubility_T
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
Clapeyron.diffusive_stability
Clapeyron.chemical_stability
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

## General Flash

```@docs
Clapeyron.FlashResult
Clapeyron.FlashData
Clapeyron.FlashSpecifications
Clapeyron.xy_flash
Clapeyron.GeneralizedXYFlash
Clapeyron.ph_flash
Clapeyron.ps_flash
Clapeyron.qt_flash
Clapeyron.qp_flash
Clapeyron.ts_flash
Clapeyron.vt_flash
```
