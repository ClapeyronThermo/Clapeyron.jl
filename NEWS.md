# v0.6.17

## New Features

- `Tproperty` and `Pproperty` speed improvements for multicomponent models.
- New method: `edge_pressure` and `edge_temperature`, that solves the isogibbs criteria for single and multicomponent models. Those functions are equivalent to `saturation_pressure`/`saturation_temperature` for single component models.
- New method: `mechanical_critical_point`, that calculates the mechanical stability limit for single and multicomponent models. For single component models, this is equivalent to `crit_pure`.
- New method: `spinodal_maximum`, that returns the maximum temperature and pressure of the diffusive spinodal line ($det(∂₂G) = 0$)

- `x0_crit_pure` now accepts an optional mol amount composition input.
- Misc documentation improvements.

## Bug Fixes

- JutulDarcy extension: fixes to allow Clapeyron work with the latest JutulDarcy extension
- CoolProp extension: fixes in `CoolProp.PropsSI` with Clapeyron models.
- CoolProp extension: fixes to support JSON parsing with CoolProp v7.
- Fixes to bubble/dew initial points.
- `SingleFluid`: Fixes when using Double exponential terms.