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
Clapeyron.spinodal_maximum
Clapeyron.edge_pressure
Clapeyron.edge_temperature
Clapeyron.UCEP_mix
Clapeyron.UCST_pressure
Clapeyron.UCST_temperature
Clapeyron.UCST_mix
Clapeyron.krichevskii_parameter
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
Clapeyron.ChemPotDewPressure
Clapeyron.FugDewPressure
Clapeyron.ActivityDewPressure
Clapeyron.ChemPotDewTemperature
Clapeyron.FugDewTemperature
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

## Flash methods

### Quick summary of available flash methods

In the following table:

- P,T represent pressure and temperature.
- X,Y represent any arbitrary specification (enthalpy, volume, internal energy, etc).
- Q represents the (molar) vapour fraction. Q = 1 is the specification of the dew point, Q = 0 corresponds to the bubble point.

| Method Name | Flash Type | maximum phases | Derivative Order | Activity Models | Electrolyte Models | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `RRTPFlash` | P-T | 2 (VLE,LLE)| First order | No (cubic EoS) | Yes | Rachford-Rice (RR) Flash |
| `MichelsenTPFlash` | P-T | 2 (VLE,LLE) | Second order¹ | Yes  | No | RR + Gibbs optimization |
| `MultiPhaseTPFlash` | P-T | any (automatic phase detection)  | Second order | Yes  | No  | RR + Gibbs optimization + phase search |
| `DETPFlash` | P-T | any (via `numphases`)  | Zeroth Order | Yes  | No | Global Gibbs optimization over fixed phases |
| `GeneralizedXYFlash` | X-Y  | 2² (VLE/LLE) | Second Order | No | No | Equilibria conditions as constrained non-linear system |
| `RRQXFlash` | Q-T,Q-P  | 2 | First Order | Yes | No | RR update + Root-finding of Q(p)/Q(T)  |
| `RRXYFlash` | P-Y,T-Y  | 2 (VLE/LLE) | First Order | Yes | No | RR iterations + Root-finding on Y(T)/Y(p) |
| `MCFlashJL` | P-T | 2 (VLE) | Second order | No | No | MultiComponentFlash.jl extension |

**Notes:**

1. `MichelsenTPFlash` can switch how it can solve the gibbs optimization problem via the `second_order::Bool` keyword, from L-BFGS (first order method) to Newton (second order Newton)  By default, `MichelsenTPFlash` uses a first order gibbs optimization solver.

1. `GeneralizedXYFlash` is a wrapper to the `xy_flash` method. While `xy_flash` can solve an arbitrary, multiphase, X-Y flash problem, `Clapeyron.jl is only capable at the moment of generating 2-phase initial points for some combinations of specifications.

The default flash methods are the following:

- Helmholtz models: `MichelsenTPFlash` (P-T), `GeneralizedXYFlash` (X-Y)
- Activity Models/Composite models: `MichelsenTPFlash` (P-T), `RRQXFlash` (Q-T,Q-P), `RRXYFlash` (P-Y,T-Y)
- Electrolyte Models: `RRTPFlash` (P-T)

### Flash API

```@docs
Clapeyron.numphases
Clapeyron.supports_reduction
Clapeyron.FlashResult
Clapeyron.FlashData
Clapeyron.FlashSpecifications
```

## Flash functions

!!! note
    For legacy reasons, the return type of `tp_flash` does not return a `FlashResult`.
    If you want to get a consistent return type for P-T flash, call `Clapeyron.tp_flash2` instead.

```@docs
Clapeyron.tp_flash
Clapeyron.ph_flash
Clapeyron.ps_flash
Clapeyron.uv_flash
Clapeyron.qt_flash
Clapeyron.qp_flash
Clapeyron.ts_flash
Clapeyron.vt_flash
Clapeyron.xy_flash
```

## General methods

```@docs
Clapeyron.DETPFlash
Clapeyron.RRTPFlash
Clapeyron.MichelsenTPFlash
Clapeyron.MultiPhaseTPFlash
Clapeyron.MCFlashJL
Clapeyron.GeneralizedXYFlash
Clapeyron.RRQXFlash
Clapeyron.RRXYFlash
```
