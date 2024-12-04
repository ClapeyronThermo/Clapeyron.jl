# v0.6.6

## New Features

- New general formulation for flashes. the formulation supports any combination of P,T,H,U,S,V,vapour fraction (q),given that initial values are provided. the formulation can be accessed via calling the function `Clapeyron.xy_flash(model,spec::FlashSpecifications,z,components0,fractions0,volumes0,T0)`
- New methods: flashes based on the general formulation, with automatic initialization:
  - p-H flash: (`ph_flash`)
  - p-S flash: (`ps_flash`)
  - V-T flash: (`vt_flash`)
  - T-S flash: (`ts_flash`)
  - vapour fraction - T flash: (`qt_flash`)
  - vapour fraction - P flash: (`qp_flash`)
- New Flash method: `GeneralizedXYFlash`, the only available method for other flashes that are not P-T formulations.
- flashes (with the exception of `tp_flash`) now return a `FlashResult` object. `Clapeyron.tp_flash2` returns a `FlashResult` that is converted to the old format internally.
- New function: `PProperty(model,T,prop,z,property)`, that calculates the pressure in T-X coordinates.
- Better `Base.show` methods for some Clapeyron.jl structs

## Bug Fixes
- `TProperty` fixes and stability improvements.
- stability improvements in calculation of bubble/dew initial points
- stability improvements when calculating Rachford-Rice iterations.