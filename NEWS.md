# v0.6.5

## New Features

- Experimental: Bulk properties for Pressure-Enthalpy and Pressure-Entropy, the syntax is the following:
  ```julia
  using Clapeyron: PH
  PH.entropy(model,p,h,z)
  PH.adiabatic_index(model,p,h,z,T0 = T0) #suplying an initial point for the temperature
  ```
  The calculation is done via `Clapeyron.Tproperty`. there are also `PT` and `VT` functions for parity.

## Bug fixes
- fixes in calculation of spinodal with cubics.
- `MultiFluid` and `SingleFluid` errors when T_reducing != Tc.
- fix `VT_identify_phase`.
