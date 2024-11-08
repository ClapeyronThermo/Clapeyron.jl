# v0.6.4

## New Features

- New model: SAFT-VR-Mie with Gross-Vrabec quatrupolar contribution (`SAFTVRMieGV`)
- New model: Co-Oriented Fluid Functional Equation for Electrostatic interactions (`COFFEE`)
- Better support for evaluation of model properties at V == Inf (ideal gas limit)
- New method: `adiabatic_index`, that calculates the ratio between the isobaric and isochoric heat capacities.
- new API: `has_fast_crit_pure`, to indicate that models can calculate their pure critical point quickly. saturation initial guesses use the result of this function to decide if and when to call the `crit_pure` routine.
- speed ups in some pressure routines
-
## Bug fixes

- `MultiFluid` and `SingleFluid` models did not use the correct gas constant.
- Fix mixing rule in `SAFTVRMie`.
- `VT_identify_phase` now returns `:unknown` for an unstable state input.
- Typos in `TProperty` for pure models.