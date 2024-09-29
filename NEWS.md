# v0.6.3

## New Features

- New Activity model: Hard-constraint Neural Network for Consistent Activity Coefficient Prediction(`HANNA`). requires the loading of the auxiliary package `ClapeyronHANNA`.
- New ideal model: PPDS correlations (`PPDSIdeal`)
- New Function: `reference_chemical_potential`, to calculate the reference chemical potential used in activity coefficient calculations
- New Functions: `activity` and `aqueous_activity`. both functions use the ability to change reference chemical potentials.
- New Functions: `spinodal_pressure` and `spinodal_temperature`.
- New function: `split_model_binaries`, that returns a list of all binary combinations of an n-component model.
- New Saturation Method: `CritExtrapolation`, that calculates saturation volumes and pressures via extrapolation from the critical point.
- `lb_volume` now has a three-arg version: `lb_volume(model,T,z)`.
- `p_scale` is now defined in terms of `T_scale` and `lb_volume`.
- Stability improvements for bubble/dew calculations.
- Speed improvements for single and multicomponent equilibria.
- Association solver is now faster for small association matrices.
- New association mixing rule: Mie-15 (`:mie15`,`:dufal`)
- Michelsen TP-Flash: in case of valid K values but single phase rachford-rice, the procedure will assume bubble or dew point as a first iteration.
- Joback: new submodule: `JobackGC` that provides all available properties using the joback correlations.
- SAFT-VR-Mie: speed improvements for calculation of association strengths.
- Cubics: Better initial point for single component saturation calculations.
- `split_model` now works for `ClapeyronParam`,`Symbol`,`Number`,`AbstractString`,`Tuple`,`Missing` and `Nothing`. before those could only be splitted if inside an `EoSModel`.
- `StructGroupParam` is deprecated, `GroupParam` has all the functionality of `StructGroupParam`.

## Bug fixes
- `SAFTgammaMie` fixes.
- `SingleFluid` has improved initial points for liquid volume evaluation.
- miscelaneous database improvements.
- `second_virial_coefficient` for cubics was ignoring the translation.
- improvements to the initial point of `SingleFluid`.